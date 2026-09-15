#!/usr/bin/env python3
"""Offline x-directed shock/turbulence diagnostics; no solver state is modified."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from nse_turbulence_statistics import _physical_local_data, _nse_indices
from slf_to_paraview_merged_cropghost import (
    read_slf, parse_step_rank, get_global_grid, get_origin_spacing,
    build_rank_range_lookup, load_meta, discover_files, group_by_step,
    complete_step_groups, select_step_groups,
)


def plane_statistics(files, meta, gamma, include_pressure=False):
    """Two-pass plane reduction; memory scales with a local SLF block, not five global fields."""
    if not np.isfinite(gamma) or gamma <= 1:
        raise ValueError('gamma must be finite and greater than one')
    first = read_slf(files[0])
    shape = tuple(get_global_grid(meta, first))
    origin, spacing = get_origin_spacing(meta, first)
    if min(shape) < 1 or shape[0] < 3 or not np.all(np.isfinite([*origin, *spacing])) or min(spacing) <= 0:
        raise ValueError('Require uniform Cartesian geometry, nx>=3 and positive spacing')
    step, _ = parse_step_rank(files[0], first)
    time = float(first.time)
    if not np.isfinite(time):
        raise ValueError('Nonfinite SLF time')
    lookup = build_rank_range_lookup(meta)
    boxes = []
    count = np.zeros(shape[0])
    sums = np.zeros((shape[0], 8))  # rho,p,u,v,w,rhou,rhov,rhow
    pressure_field = np.empty(shape) if include_pressure else None

    def blocks():
        for path in files:
            slf = read_slf(path)
            s, rank = parse_step_rank(path, slf)
            if s != step or not np.isclose(slf.time, time, rtol=1e-12, atol=1e-14):
                raise ValueError('Mixed step or time in rank files')
            if tuple(int(n) for n in slf.meta[2:5]) != shape:
                raise ValueError('SLF global grid and metadata disagree')
            if rank is not None and not lookup:
                raise ValueError('Rank layout requires explicit rank_ranges in meta.json')
            data, box = _physical_local_data(slf, rank, lookup, shape)
            if any(s.start < 0 or s.stop > n or s.stop <= s.start for s, n in zip(box, shape)):
                raise ValueError('Rank range outside global grid')
            if len(slf.names) != 5 or not np.all(np.isfinite(data)):
                raise ValueError('Require finite single-component NSE state (five variables)')
            ri, ui, vi, wi, ei = _nse_indices(slf.names)
            rho = data[..., ri]
            if np.any(rho <= 0):
                raise ValueError('Nonpositive density; no clipping')
            vel = data[..., [ui, vi, wi]] / rho[..., None]
            p = (gamma-1)*(data[..., ei] - .5*rho*np.sum(vel**2, axis=-1))
            if not np.all(np.isfinite(p)) or np.any(p <= 0):
                raise ValueError('Invalid pressure; no clipping')
            yield box, rho, p, vel

    for box, rho, p, vel in blocks():
        if any(all(a.start < b.stop and b.start < a.stop for a, b in zip(box, other)) for other in boxes):
            raise ValueError('Overlapping rank blocks')
        boxes.append(box)
        if include_pressure:
            pressure_field[box] = p
        ix = box[0]
        count[ix] += rho.shape[1]*rho.shape[2]
        sums[ix, 0] += rho.sum(axis=(1, 2))
        sums[ix, 1] += p.sum(axis=(1, 2))
        sums[ix, 2:5] += vel.sum(axis=(1, 2))
        sums[ix, 5:8] += (rho[..., None]*vel).sum(axis=(1, 2))
    if not np.all(count == shape[1]*shape[2]):
        raise ValueError('Missing cells in rank layout')
    means = sums[:, 2:5]/count[:, None]
    favre = sums[:, 5:8]/sums[:, 0, None]
    cov = np.zeros((shape[0], 6)); fcov = np.zeros_like(cov)
    pairs = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
    for box, rho, p, vel in blocks():
        ix = box[0]
        du = vel-means[ix, None, None, :]
        df = vel-favre[ix, None, None, :]
        for j, (a, b) in enumerate(pairs):
            cov[ix, j] += (du[..., a]*du[..., b]).sum(axis=(1,2))
            fcov[ix, j] += (rho*df[..., a]*df[..., b]).sum(axis=(1,2))
    cov /= count[:, None]; fcov /= sums[:, 0, None]
    profile = dict(x=origin[0]+(np.arange(shape[0])+.5)*spacing[0],
                   rho=sums[:,0]/count, p=sums[:,1]/count)
    for j, name in enumerate('uvw'):
        profile['mean_'+name] = means[:,j]
        profile['favre_'+name] = favre[:,j]
    for j, name in enumerate(('uu','vv','ww','uv','uw','vw')):
        profile['R_'+name] = cov[:,j]
        profile['F_'+name] = fcov[:,j]
    profile['tke'] = .5*cov[:,:3].sum(axis=1)
    profile['favre_tke'] = .5*fcov[:,:3].sum(axis=1)
    result = (step, time, profile, spacing)
    return (*result, pressure_field) if include_pressure else result


def shock_surface(pressure, x, direction, search, upstream, downstream,
                  min_ratio=1.05, previous=None, max_shift=None):
    """Independent x-line detection. Invalid lines return NaN, never the mean front."""
    pressure = np.asarray(pressure)
    x = np.asarray(x)
    if pressure.ndim != 3 or len(x) != pressure.shape[0] or len(x) < 3:
        raise ValueError('Invalid pressure field/x axis')
    if not np.all(np.isfinite(pressure)) or np.any(pressure <= 0) or np.any(np.diff(x) <= 0):
        raise ValueError('Require positive finite pressure and increasing x')
    shape = pressure.shape[1:]
    if previous is not None and np.shape(previous) != shape:
        raise ValueError('Previous shock surface grid mismatch')
    position = np.full(shape, np.nan)
    strength = np.zeros(shape)
    # Only 2D work arrays: avoid constructing a second global pressure-sized gradient field.
    for i, face in enumerate((x[:-1]+x[1:])/2):
        if not search[0] <= face <= search[1]:
            continue
        gradient = -direction*(pressure[i+1]-pressure[i])/(x[i+1]-x[i])
        choose = gradient > strength
        if previous is not None and max_shift is not None:
            choose &= ~np.isfinite(previous) | (abs(face-previous) <= max_shift)
        position[choose] = face
        strength[choose] = gradient[choose]
    valid = np.isfinite(position)
    averages = []
    dx = x[1]-x[0]
    for limits, sign in ((upstream,1),(downstream,-1)):
        count = np.zeros(shape, dtype=np.int64); total = np.zeros(shape)
        for i, coordinate in enumerate(x):
            distance = sign*direction*(coordinate-position)
            inside = (distance >= limits[0]) & (distance <= limits[1])
            count += inside; total += np.where(inside, pressure[i], 0)
        extent = np.maximum(sign*direction*(x[0]-dx/2-position),
                            sign*direction*(x[-1]+dx/2-position))
        valid &= (count >= 2) & (extent >= limits[1]-1e-10*dx)
        averages.append(np.divide(total,count,out=np.full(shape,np.nan),where=count>0))
    ratio = averages[1]/averages[0]
    valid &= np.isfinite(ratio) & (ratio >= min_ratio)
    position[~valid] = np.nan
    strength[~valid] = np.nan
    return dict(shock_x=position, valid=valid, pressure_gradient=strength, pressure_ratio=ratio)


def surface_summary(surface):
    values = surface['shock_x'][surface['valid']]
    result = dict(surface_valid_count=int(values.size), surface_total_count=int(surface['valid'].size),
                  surface_valid_fraction=float(np.mean(surface['valid'])))
    for key, function in [('mean',np.mean),('std',np.std),('min',np.min),('max',np.max)]:
        result['surface_x_'+key] = float(function(values)) if values.size else float('nan')
    return result


def diagnose(profile, direction, search, upstream, downstream, previous=None, max_shift=None, min_ratio=1.05):
    x, p = profile['x'], profile['p']
    face = (x[1:]+x[:-1])/2
    gradient = -direction*np.diff(p)/np.diff(x)
    allowed = (face >= search[0]) & (face <= search[1]) & (gradient > 0)
    if previous is not None and max_shift is not None:
        allowed &= abs(face-previous) <= max_shift
    if not allowed.any():
        raise ValueError('No compressive pressure front within search/tracking range')
    index = int(np.argmax(np.where(allowed, gradient, -np.inf)))
    shock = float(face[index])
    distance = direction*(x-shock)
    result = dict(shock_x=shock, pressure_gradient=float(gradient[index]))
    for name, limits, sign in [('upstream',upstream,1),('downstream',downstream,-1)]:
        signed = sign*distance
        # Reject truncated windows, rather than silently changing sample regions near boundaries.
        dx = x[1]-x[0]
        edges = sign*direction*(np.array([x[0]-dx/2,x[-1]+dx/2])-shock)
        if limits[1] > edges.max()+1e-10*dx:
            raise ValueError(name+' window extends outside domain')
        mask = (signed >= limits[0]) & (signed <= limits[1])
        if mask.sum() < 2:
            raise ValueError(name+' window needs at least two planes')
        result[name+'_planes'] = int(mask.sum())
        for key, values in profile.items():
            if key == 'x':
                continue
            weights = profile['rho'][mask] if key.startswith(('favre_', 'F_')) else None
            result[name+'_'+key] = float(np.average(values[mask], weights=weights))
    ratio = result['downstream_p']/result['upstream_p']
    if ratio < min_ratio:
        raise ValueError('Pressure ratio below threshold; check front identity and sample windows')
    result['pressure_ratio'] = ratio
    for key in ('tke','favre_tke','R_uu','R_vv','R_ww','F_uu','F_vv','F_ww'):
        denominator = result['upstream_'+key]
        result['amplification_'+key] = result['downstream_'+key]/denominator if denominator > 1e-30 else float('nan')
    return result


def write_csv(path, rows):
    with path.open('x', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--meta', type=Path, required=True)
    parser.add_argument('--gamma', type=float, required=True)
    parser.add_argument('--direction', type=int, choices=(-1,1), required=True)
    parser.add_argument('--search', type=float, nargs=2, required=True, metavar=('XMIN','XMAX'))
    parser.add_argument('--upstream', type=float, nargs=2, required=True, metavar=('NEAR','FAR'))
    parser.add_argument('--downstream', type=float, nargs=2, required=True, metavar=('NEAR','FAR'))
    parser.add_argument('--max-shift', type=float, help='Maximum front displacement between selected outputs')
    parser.add_argument('--min-pressure-ratio', type=float, default=1.05)
    parser.add_argument('--steps', default='all')
    parser.add_argument('--layout', choices=('auto','global','rank'), default='auto')
    args = parser.parse_args(argv)
    try:
        for limits in (args.search, args.upstream, args.downstream):
            if not np.all(np.isfinite(limits)) or limits[1] <= limits[0]:
                raise ValueError('Window limits must be finite and increasing')
        if min(args.upstream[0],args.downstream[0]) <= 0:
            raise ValueError('NEAR distances must be positive to exclude the shock')
        if not np.isfinite(args.min_pressure_ratio) or args.min_pressure_ratio <= 1:
            raise ValueError('min-pressure-ratio must exceed one')
        if args.max_shift is not None and (not np.isfinite(args.max_shift) or args.max_shift <= 0):
            raise ValueError('max-shift must be positive')
        meta = load_meta(args.meta)
        if str(meta.get('equation','nse')).lower() != 'nse':
            raise ValueError('Only single-component NSE is supported')
        files = discover_files(args.input, case_name=meta.get('case_name') or None, layout=args.layout, meta=meta)
        groups = select_step_groups(complete_step_groups(group_by_step(files),meta),args.steps)
        if not groups:
            raise ValueError('No complete selected SLF steps')
        args.output_dir.mkdir(parents=True, exist_ok=False)
        records = []; previous = None; geometry = None; previous_surface = None
        for step in sorted(groups):
            step_read,time,profile,spacing,pressure = plane_statistics(groups[step],meta,args.gamma,include_pressure=True)
            if records and time <= records[-1]['time']:
                raise ValueError('Selected times must be strictly increasing')
            if geometry is not None and not np.array_equal(profile['x'],geometry):
                raise ValueError('Changing grids are unsupported')
            geometry = profile['x'].copy()
            record = dict(step=step_read,time=time,**diagnose(profile,args.direction,args.search,
                args.upstream,args.downstream,previous,args.max_shift,args.min_pressure_ratio))
            surface = shock_surface(pressure,profile['x'],args.direction,args.search,args.upstream,args.downstream,
                                    args.min_pressure_ratio,previous_surface,args.max_shift)
            del pressure
            record.update(surface_summary(surface))
            print(f"[INFO] step={step_read} surface valid={record['surface_valid_count']}/{record['surface_total_count']}")
            previous_surface = surface['shock_x'].copy()
            first = read_slf(groups[step][0])
            origin,_ = get_origin_spacing(meta,first)
            y = origin[1]+(np.arange(surface['shock_x'].shape[0])+.5)*spacing[1]
            z = origin[2]+(np.arange(surface['shock_x'].shape[1])+.5)*spacing[2]
            del first
            np.savez_compressed(args.output_dir/f'shock_surface_{step_read:08d}.npz',
                                step=step_read,time=time,y=y,z=z,**surface)
            previous = record['shock_x']; records.append(record)
            write_csv(args.output_dir/f'planes_{step_read:08d}.csv',
                [dict(zip(profile,values)) for values in zip(*profile.values())])
        times = np.array([r['time'] for r in records]); positions = np.array([r['shock_x'] for r in records])
        speeds = np.gradient(positions,times) if len(records)>1 else [float('nan')]
        for row,speed in zip(records,speeds):
            row['shock_speed'] = float(speed)
        write_csv(args.output_dir/'shock_history.csv',records)
        definitions = dict(status='complete', parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
            source_files=[str(p.resolve()) for s in sorted(groups) for p in groups[s]],
            definitions={
                'shock_x':'Face of largest signed gradient of yz-mean pressure within search/tracking range; grid-scale resolution.',
                'shock_surface':'NPZ shock_x[j,k]=x_s(y[j],z[k],time), independent pressure-gradient maximum on each x line; invalid lines are NaN with valid=False.',
                'surface_statistics':'Mean, population standard deviation (ddof=0), min/max over valid lines; equal transverse cell areas. Always inspect valid_fraction.',
                'surface_tracking':'max-shift applied per line to previous selected valid position; invalid previous lines reacquire within global search range.',
                'shock_speed':'Signed dx_shock/dt using saved times; interior nonuniform central differences, endpoints one-sided. Single time: NaN.',
                'R':'Plane Reynolds covariance about each yz-plane volume mean; uu,vv,ww,uv,uw,vw.',
                'F':'Plane density-weighted covariance about each yz-plane Favre mean.',
                'tke':'0.5*(R_uu+R_vv+R_ww), per unit mass; Favre analog uses F.',
                'regions':'Distances from tracked front, upstream in propagation direction. Average plane statistics, not whole-region velocity variance.',
                'amplification':'Simultaneous downstream/upstream ratio, not matched fluid parcels. Denominator <=1e-30 yields NaN.',
                'units':'Same coordinates, time and velocity normalization as SLF; no SI conversion.',
                'limitations':'Uniform Cartesian x-directed single-valued front only. Mean-profile and linewise fronts are separate. Region turbulence statistics still use the mean-profile front, not a surface-fitted window. Search must isolate shock from other waves.'})
        (args.output_dir/'analysis.json').write_text(json.dumps(definitions,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print('[OK] Shock analysis completed:',args.output_dir)
        return 0
    except (OSError,ValueError,EOFError) as exc:
        print('[ERROR]',exc)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
