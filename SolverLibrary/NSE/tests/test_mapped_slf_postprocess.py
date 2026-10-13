"""Mapped SLF geometry, rank integrity, binary VTK and weighted-statistics tests."""
import math
import struct
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import slf_to_paraview_merged_cropghost as pv
import nse_turbulence_statistics as stats


def fixture(root, rank=0, count=1, norm=1, periodic=1):
    n, ng = 6, 3
    js, je = (1, 6) if count == 1 else (1+3*rank, 3+3*rank)
    h = (1, 10, n, n, n, js, je, 1, n, norm, 2 if norm else 1, 1, periodic, periodic, periodic)
    r = (.25, 1.67, 100., .72, 1., .1, 0., 2., -1., 1., 0., 3., .7, 1., 1.5)
    shape = (n+2*ng, je-js+1+2*ng, n+2*ng, 5)
    meta = (10, rank, n, n, n, ng, count, 1)
    dummy = pv.SLFData(Path(''), 1, 2, shape, np.array(meta), r[0], r[6:12], [], None, h, r)
    edges, _, _ = pv.mapped_axes(dummy)
    x = (edges[0][1:]+edges[0][:-1])/2
    q = np.full(shape, -999., dtype='<f8')
    interior = q[ng:-ng, ng:-ng, ng:-ng]
    interior[..., 0] = 2
    interior[..., 1] = 2*x[:, None, None]
    interior[..., 2:4] = 0
    interior[..., 4] = 3/(r[1]-1)+x[:, None, None]**2
    path = root/f'field_000010_rank{rank:05d}.slf'
    with path.open('wb') as f:
        f.write(b'SLF1\0\0\0\0'+struct.pack('<3i4i8id6di', 1, 2, 4, *shape, *meta, r[0], *r[6:12], 5))
        for name in ('rho', 'rho_u', 'rho_v', 'rho_w', 'rho_E'):
            f.write(name.encode().ljust(32, b' '))
        f.write(q.tobytes(order='F'))
        f.write(b'NSEPAR1\0'+struct.pack('<5d', *r[1:6])+('central6' if norm else 'fv2').encode().ljust(32, b' '))
        f.write(b'NSEGRID1'+struct.pack('<15i15d', *h, *r)+b'NSREND1\0')
    return path


class MappedSlfTests(unittest.TestCase):
    def test_rank_merge_and_weighted_statistics(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [fixture(Path(tmp), rank, 2) for rank in range(2)]
            first, q = pv.merge_mapped(paths)
            self.assertEqual(q.shape, (6, 6, 6, 5))
            edges, volume, weight = pv.mapped_axes(first)
            self.assertAlmostEqual(volume.sum(), 12.)
            self.assertFalse(np.allclose(volume, weight))
            record = stats.mapped_statistics(paths, 1e-12, 1e-12)
            self.assertAlmostEqual(record['mean_density'], 2.)
            self.assertAlmostEqual(record['mean_pressure'], 3.)
            self.assertAlmostEqual(record['mean_u'], 1.)
            self.assertAlmostEqual(record['volume_integral_mass'], 24.)
            self.assertAlmostEqual(record['solver_integral_mass'], 2*weight.sum())
            expected = np.sum(((edges[0][1:]+edges[0][:-1])/2-1)**2*np.diff(edges[0]))/2
            self.assertAlmostEqual(record['rms_u_fluctuation']**2, expected)
            self.assertTrue(math.isnan(record['dissipation_rate']))

    def test_vtr_binary_coordinates_values_and_pvd(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            self.assertEqual(pv.main([str(root), '-o', str(root/'pv'), '--fields', 'p,u']), 0)
            raw = (root/'pv/field_000010.vtr').read_bytes()
            prefix, tail = raw.split(b'<AppendedData encoding="raw">_', 1)
            xml = ET.fromstring(prefix+b'</VTKFile>')
            arrays = xml.findall('.//CellData/DataArray')+xml.findall('.//Coordinates/DataArray')
            values = {}
            for a in arrays:
                offset = int(a.attrib['offset'])
                size = struct.unpack_from('<Q', tail, offset)[0]
                values[a.attrib['Name']] = np.frombuffer(tail[offset+8:offset+8+size], dtype='<f8')
            np.testing.assert_allclose(values['p'], 3.)
            first = pv.read_slf(root/'field_000010_rank00000.slf', allow_nonuniform=True)
            edges, volume, weight = pv.mapped_axes(first)
            for a, name in enumerate('xyz'):
                np.testing.assert_allclose(values[name], edges[a])
            np.testing.assert_allclose(values['integration_weight'], weight.ravel(order='F'))
            pvd = ET.parse(root/'pv/collection.pvd')
            self.assertEqual(pvd.find('.//DataSet').attrib['timestep'], '0.25')

    def test_missing_duplicate_and_mixed_ranks_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, b = fixture(root, 0, 2), fixture(root, 1, 2)
            for files in ([a], [a, a]):
                with self.assertRaises(ValueError):
                    pv.merge_mapped(files)
            fixture(root, 1, 2, periodic=0)
            with self.assertRaises(ValueError):
                pv.merge_mapped([a, b])

    def test_truncation_and_fft_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = fixture(Path(tmp))
            with self.assertRaises(ValueError):
                pv.read_slf(p)
            p.write_bytes(p.read_bytes()[:-1])
            with self.assertRaises((ValueError, EOFError)):
                pv.read_slf(p, allow_nonuniform=True)

    def test_fv2_weights_and_statistics_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = fixture(root, norm=0, periodic=0)
            first = pv.read_slf(p, allow_nonuniform=True)
            _, volume, weight = pv.mapped_axes(first)
            np.testing.assert_array_equal(volume, weight)
            self.assertEqual(stats.main([str(root), '-o', str(root/'stats.csv'), '--gamma', '1.4', '--reynolds', '999']), 0)
            self.assertIn('solver_integral_mass', (root/'stats.csv').read_text())
            self.assertIn('no FFT', (root/'stats_metadata.json').read_text())

    def test_stride_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture(root)
            with self.assertRaises(SystemExit):
                pv.main([str(root), '-o', str(root/'pv'), '--stride', '2'])
            self.assertFalse((root/'pv').exists())


if __name__ == '__main__':
    unittest.main()
