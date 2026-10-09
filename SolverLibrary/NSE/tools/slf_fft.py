#!/usr/bin/env python3
"""Spatial FFT of a scalar SLF field; rank assembly is independent of MPI size."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from slf_to_paraview_merged_cropghost import (
    read_slf, load_meta, get_global_grid, get_origin_spacing,
    build_rank_range_lookup, rank_slices, crop_center_to_shape,
    discover_files, group_by_step, select_step_groups,
)


def scalar_field(slf, field, gamma):
    lookup = {name.lower(): i for i, name in enumerate(slf.names)}
    if field in lookup:
        return slf.data[..., lookup[field]]
    if field not in ('u', 'v', 'w', 'p'):
        raise ValueError(f"Unknown field {field!r}; stored fields: {slf.names}")
    rho = slf.data[..., lookup['rho']]
    if np.any(rho <= 0):
        raise ValueError(f"{slf.path}: nonpositive density")
    if field != 'p':
        return slf.data[..., lookup['rho_' + field]] / rho
    if gamma is None or not np.isfinite(gamma) or gamma <= 1:
        raise ValueError('Derived p requires --gamma (>1), matching the single-component ideal gas case')
    kinetic = sum(slf.data[..., lookup['rho_' + c]] ** 2 for c in 'uvw') / (2 * rho)
    return (gamma - 1) * (slf.data[..., lookup['rho_e']] - kinetic)


def assemble(files, meta, field, gamma=None):
    """Reject incomplete/overlapping blocks rather than silently FFT bad data."""
    first = read_slf(files[0])
    shape = get_global_grid(meta, first)
    header_shape = tuple(int(n) for n in first.meta[2:5])
    if all(n > 0 for n in header_shape) and tuple(shape) != header_shape:
        raise ValueError('Global grid differs from SLF header; supply matching meta.json')
    if any(n <= 0 for n in shape):
        raise ValueError('Invalid global grid')
    lookup = build_rank_range_lookup(meta)
    result = np.empty(shape)
    covered = np.zeros(shape, dtype=bool)
    seen = set()
    for index, path in enumerate(files):
        slf = first if index == 0 else read_slf(path)
        rank = int(slf.meta[1])
        if rank in seen:
            raise ValueError('Duplicate rank / mixed global and rank files')
        seen.add(rank)
        if slf.time != first.time or slf.meta[0] != first.meta[0] or slf.names != first.names:
            raise ValueError('Inconsistent step, time or fields across ranks')
        if slf.bounds != first.bounds:
            raise ValueError('Inconsistent domain bounds across ranks')
        distributed = int(slf.meta[6]) > 1 or len(files) > 1
        if distributed and rank not in lookup:
            raise ValueError('Distributed SLF requires matching meta.json rank_ranges')
        slices = rank_slices(rank, slf, lookup) if rank in lookup else tuple(slice(0, n) for n in shape)
        if any(s.start < 0 or s.stop > n or s.stop <= s.start for s, n in zip(slices, shape)):
            raise ValueError('Invalid rank bounds')
        if covered[slices].any():
            raise ValueError('Overlapping rank blocks')
        expected = tuple(s.stop - s.start for s in slices)
        # Derive primitives after removing ghosts (which may be uninitialized).
        slf.data = crop_center_to_shape(slf.data, expected, path)
        local = scalar_field(slf, field, gamma)
        if not np.isfinite(local).all():
            raise ValueError(f'{path}: NaN/Inf in requested field')
        result[slices] = local
        covered[slices] = True
    if not covered.all():
        raise ValueError('Missing rank blocks: global field is incomplete')
    _, spacing = get_origin_spacing(meta, first)
    if not np.isfinite(spacing).all() or min(spacing) <= 0:
        raise ValueError('Invalid spacing')
    return result, spacing, float(first.time)


def analyze(data, spacing, remove_mean=True, window='none'):
    """F=fftn(f)/N; sum |F|^2 = mean(f^2), without one-sided factors."""
    field = np.array(data, dtype=float, copy=True)
    mean = float(field.mean())
    if remove_mean:
        field -= mean
    if window == 'hann':
        for axis, n in enumerate(field.shape):
            if n < 3:
                raise ValueError('Hann window requires at least three cells per axis')
            shape = [1, 1, 1]
            shape[axis] = n
            field *= np.hanning(n).reshape(shape)
    elif window != 'none':
        raise ValueError('Unknown window')
    transform = np.fft.fftn(field) / field.size
    power = np.abs(transform) ** 2
    axes = [2 * np.pi * np.fft.fftfreq(n, d=h) for n, h in zip(field.shape, spacing)]
    # Physical angular wave numbers; rectangular domains are supported.
    dk = min(2 * np.pi / (n * h) for n, h in zip(field.shape, spacing))
    radius = np.sqrt(axes[0][:, None, None]**2 + axes[1][None, :, None]**2 + axes[2][None, None, :]**2)
    bins = np.floor(radius / dk + 0.5).astype(np.int64)
    counts = np.bincount(bins.ravel())
    shell = np.bincount(bins.ravel(), weights=power.ravel(), minlength=len(counts))
    table = np.column_stack((np.arange(len(counts))*dk, counts, shell, shell/dk))
    info = dict(mean=mean, mean_removed=remove_mean, window=window,
                physical_mean_square=float(np.mean(field**2)), spectral_sum=float(power.sum()),
                shell_width=dk, normalization='F=fftn(processed_field)/N; sum(abs(F)^2)=mean(processed_field^2)',
                window_power_corrected=False)
    return transform, axes, table, info


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='SLF file or output directory')
    parser.add_argument('-o', '--output-dir', type=Path, required=True)
    parser.add_argument('--meta', type=Path)
    parser.add_argument('--field', default='rho', help='stored field or NSE u,v,w,p (one field per invocation)')
    parser.add_argument('--gamma', type=float, help='required when deriving pressure')
    parser.add_argument('--steps', default='latest', help='latest, all, 100, or 0:1000:100')
    parser.add_argument('--layout', choices=('auto', 'global', 'rank'), default='auto')
    parser.add_argument('--keep-mean', action='store_true')
    parser.add_argument('--window', choices=('none', 'hann'), default='none')
    parser.add_argument('--save-fft', action='store_true', help='also save full complex FFT and angular wave-number axes as NPZ')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args(argv)
    try:
        field = args.field.lower()
        if not field.replace('_', '').isalnum():
            raise ValueError('Invalid field name')
        directory = args.input if args.input.is_dir() else args.input.parent
        meta_path = args.meta or directory / 'meta.json'
        meta = load_meta(meta_path if meta_path.exists() or args.meta else None)
        files = discover_files(args.input, layout=args.layout, meta=meta)
        groups = select_step_groups(group_by_step(files), args.steps)
        if not groups:
            raise ValueError('No SLF files found')
        print('NOTE: spatial FFT assumes a uniform Cartesian grid and periodic continuation; no temporal FFT.')
        for step, paths in groups.items():
            stem = args.output_dir / f'fft_{step:06d}_{field}'
            outputs = [stem.with_suffix('.csv'), stem.with_suffix('.json')]
            if args.save_fft:
                outputs.append(stem.with_suffix('.npz'))
            if not args.overwrite and any(p.exists() for p in outputs):
                raise ValueError(f'Output exists: {stem}; use --overwrite explicitly')
            data, spacing, time = assemble(paths, meta, field, args.gamma)
            transform, axes, table, info = analyze(data, spacing, not args.keep_mean, args.window)
            info.update(step=step, time=time, field=field, gamma=args.gamma, shape=list(data.shape),
                        spacing=list(spacing), sources=[str(p.resolve()) for p in paths])
            args.output_dir.mkdir(parents=True, exist_ok=True)
            np.savetxt(outputs[0], table, delimiter=',', comments='',
                       header='k,mode_count,shell_power,power_per_unit_k')
            outputs[1].write_text(json.dumps(info, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            if args.save_fft:
                np.savez(outputs[2], fft=transform, kx=axes[0], ky=axes[1], kz=axes[2])
            print(f'[OK] {outputs[0]} Parseval: {info["physical_mean_square"]:.12g} / {info["spectral_sum"]:.12g}')
        return 0
    except (ValueError, KeyError, OSError, EOFError) as exc:
        print(f'[ERROR] {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
