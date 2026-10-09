from pathlib import Path
import json
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from slf_fft import analyze, assemble, main


def write_slf(path, data, rank=0, nprocs=1, global_shape=None, time=0.1):
    names = ['rho', 'rho_u', 'rho_v', 'rho_w', 'rho_E'][:data.shape[-1]]
    meta = np.array([10, rank, *(global_shape or data.shape[:3]), 0, nprocs, 0], dtype='<i4')
    with path.open('wb') as f:
        f.write(b'SLF1\0\0\0\0')
        f.write(struct.pack('<iii', 1, 2, 3))
        f.write(np.array(data.shape, dtype='<i4').tobytes())
        f.write(meta.tobytes())
        f.write(struct.pack('<d6di', time, 0, 8, 0, 6, 0, 4, len(names)))
        for name in names:
            f.write(name.encode().ljust(32, b'\0'))
        f.write(np.asarray(data, dtype='<f8').tobytes(order='F'))


class FFTTests(unittest.TestCase):
    def test_sine_parseval_and_wave_number(self):
        data = np.broadcast_to(3 + np.sin(2*np.pi*2*np.arange(8)/8)[:, None, None], (8, 6, 4))
        fft, axes, table, info = analyze(data, (1, 1, 1))
        self.assertAlmostEqual(abs(fft[2, 0, 0]), 0.5)
        self.assertAlmostEqual(axes[0][2], np.pi/2)
        self.assertAlmostEqual(table[:, 2].sum(), 0.5)
        self.assertAlmostEqual(info['physical_mean_square'], info['spectral_sum'])
        self.assertAlmostEqual(info['mean'], 3)

    def test_constant_mean_and_hann(self):
        data = np.ones((8, 7, 6))
        self.assertAlmostEqual(analyze(data, (1, 2, 3))[3]['spectral_sum'], 0)
        self.assertAlmostEqual(analyze(data, (1, 2, 3), False)[3]['spectral_sum'], 1)
        info = analyze(data, (1, 2, 3), False, 'hann')[3]
        self.assertAlmostEqual(info['physical_mean_square'], info['spectral_sum'])
        self.assertLess(info['spectral_sum'], 1)

    def test_rank_merge_ghost_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = np.arange(8*6*4).reshape(8, 6, 4, 1).astype(float)
            paths = [root / f'field_000010_rank{r:05d}.slf' for r in range(2)]
            ranges = []
            for r, path in enumerate(paths):
                local = np.pad(data[:, r*3:(r+1)*3], ((1,1),(1,1),(1,1),(0,0)))
                write_slf(path, local, r, 2, (8,6,4))
                ranges.append(dict(rank=r, i_start=1, i_end=8, j_start=1+r*3, j_end=3+r*3, k_start=1, k_end=4))
            meta = dict(grid=[8,6,4], parallel=dict(rank_ranges=ranges))
            merged, spacing, _ = assemble(paths, meta, 'rho')
            np.testing.assert_array_equal(merged, data[...,0])
            self.assertEqual(spacing, (1,1,1))
            with self.assertRaises(ValueError): assemble(paths[:1], meta, 'rho')
            with self.assertRaises(ValueError): assemble(paths, {}, 'rho')
            with self.assertRaises(ValueError): assemble(paths + paths[:1], meta, 'rho')
            write_slf(paths[1], data[:, 3:], 1, 2, (8,6,4), time=0.2)
            with self.assertRaises(ValueError): assemble(paths, meta, 'rho')

    def test_pressure_and_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'field_000010.slf'
            data = np.zeros((8,6,4,5))
            data[...,0] = 2
            data[...,1] = 6
            data[...,4] = 14
            write_slf(path, data)
            p, _, _ = assemble([path], {}, 'p', 1.4)
            np.testing.assert_allclose(p, 2)
            with self.assertRaises(ValueError): assemble([path], {}, 'p')
            args = [str(path), '-o', str(root/'fft'), '--field', 'p', '--gamma', '1.4', '--save-fft']
            self.assertEqual(main(args), 0)
            saved = np.load(root/'fft/fft_000010_p.npz')
            self.assertEqual(saved['fft'].shape, (8,6,4))
            saved.close()
            info = json.loads((root/'fft/fft_000010_p.json').read_text())
            self.assertAlmostEqual(info['mean'], 2)
            self.assertEqual(main(args), 1)  # No silent overwrite.
            self.assertEqual(main(args + ['--overwrite']), 0)
            data[0,0,0,0] = np.nan
            write_slf(path, data)
            with self.assertRaises(ValueError): assemble([path], {}, 'rho')


if __name__ == '__main__':
    unittest.main()
