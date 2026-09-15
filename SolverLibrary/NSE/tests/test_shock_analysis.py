import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from nse_shock_analysis import plane_statistics, diagnose, main, shock_surface, surface_summary
def write_slf(path,data,rank=0,nprocs=1,global_shape=None,time=.1):
    import struct
    names=['rho','rho_u','rho_v','rho_w','rho_E']
    with path.open('wb') as f:
        f.write(b'SLF1\0\0\0\0');f.write(struct.pack('<iii',1,2,3))
        f.write(np.array(data.shape,dtype='<i4').tobytes())
        f.write(np.array([10,rank,*(global_shape or data.shape[:3]),0,nprocs,0],dtype='<i4').tobytes())
        f.write(struct.pack('<d6di',time,0,40,0,4,0,4,5))
        for name in names:f.write(name.encode().ljust(32,b'\0'))
        f.write(np.asarray(data,dtype='<f8').tobytes(order='F'))


def state(nx=40,shock=20,direction=1):
    x=(np.arange(nx)+.5)[:,None,None]
    behind=direction*(x-shock)<0
    rho=np.broadcast_to(np.where(behind,2.,1.),(nx,4,4))
    p=np.broadcast_to(np.where(behind,4.,1.),rho.shape)
    fluct=np.array([-1.,1.,-1.,1.])[None,:,None]
    # Nonuniform mean along x must NOT be counted as turbulent energy.
    u=np.broadcast_to(.1*x+np.where(behind,2.,1.)*fluct,rho.shape)
    q=np.zeros((*rho.shape,5));q[...,0]=rho;q[...,1]=rho*u
    q[...,4]=p/.4+.5*rho*u*u
    return q


class ShockTests(unittest.TestCase):
    def test_corrugated_surface_both_directions_and_missing_lines(self):
        x=np.arange(40)+.5
        exact=np.array([[17,18,19,20],[18,19,20,21],[19,20,21,22],[20,21,22,23]],dtype=float)
        for direction in (1,-1):
            pressure=np.where(direction*(x[:,None,None]-exact[None,:,:])<0,4.,1.)
            found=shock_surface(pressure,x,direction,[10,30],[2,8],[2,8])
            np.testing.assert_array_equal(found['shock_x'],exact)
            stats=surface_summary(found)
            self.assertEqual(stats['surface_x_mean'],exact.mean())
            self.assertEqual(stats['surface_x_std'],exact.std())
            self.assertEqual(stats['surface_valid_fraction'],1)
            pressure[:,0,0]=1
            missing=shock_surface(pressure,x,direction,[10,30],[2,8],[2,8])
            self.assertFalse(missing['valid'][0,0])
            self.assertTrue(np.isnan(missing['shock_x'][0,0]))
            self.assertEqual(surface_summary(missing)['surface_valid_count'],15)
            restricted=shock_surface(pressure,x,direction,[10,30],[2,8],[2,8],previous=exact-5,max_shift=1)
            self.assertFalse(restricted['valid'].any())

    def test_x_partition_surface_matches_global(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);q=state();ranges=[];paths=[]
            for rank in range(2):
                path=root/f'field_000010_rank{rank:05d}.slf';paths.append(path)
                write_slf(path,np.pad(q[rank*20:(rank+1)*20],((1,1),(1,1),(1,1),(0,0))),rank,2,(40,4,4))
                ranges.append(dict(rank=rank,i_start=rank*20+1,i_end=(rank+1)*20,j_start=1,j_end=4,k_start=1,k_end=4))
            meta=dict(grid=[40,4,4],spacing=[1,1,1],parallel=dict(rank_ranges=ranges))
            _,_,profile,_,pressure=plane_statistics(paths,meta,1.4,include_pressure=True)
            found=shock_surface(pressure,profile['x'],1,[10,30],[2,8],[2,8])
            np.testing.assert_array_equal(found['shock_x'],np.full((4,4),20.))
    def test_plane_statistics_and_amplification(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'field_000010.slf';write_slf(p,state())
            _,_,prof,_=plane_statistics([p],dict(grid=[40,4,4],spacing=[1,1,1]),1.4)
            r=diagnose(prof,1,[10,30],[2,8],[2,8])
            self.assertEqual(r['shock_x'],20)
            self.assertAlmostEqual(r['amplification_tke'],4)
            self.assertAlmostEqual(r['amplification_favre_tke'],4)
            self.assertAlmostEqual(r['upstream_tke'],.5)
            with self.assertRaises(ValueError):diagnose(prof,1,[10,30],[2,80],[2,8])
            with self.assertRaises(ValueError):diagnose(prof,-1,[10,30],[2,8],[2,8])

    def test_favre_and_rank_equivalence(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);q=state();rho=q[...,0].copy();u=q[...,1]/rho
            p=.4*(q[...,4]-.5*rho*u*u)
            q[...,0]*=np.array([1.,2.,1.,2.])[None,:,None]
            q[...,1]=q[...,0]*u;q[...,4]=p/.4+.5*q[...,0]*u*u
            whole=root/'field_000010.slf';write_slf(whole,q)
            meta=dict(grid=[40,4,4],spacing=[1,1,1])
            expected=plane_statistics([whole],meta,1.4)[2]
            paths=[];ranges=[]
            for rank in range(2):
                path=root/f'field_000010_rank{rank:05d}.slf';paths.append(path)
                write_slf(path,np.pad(q[:,rank*2:(rank+1)*2],((1,1),(1,1),(1,1),(0,0))),rank,2,(40,4,4))
                ranges.append(dict(rank=rank,i_start=1,i_end=40,j_start=rank*2+1,j_end=rank*2+2,k_start=1,k_end=4))
            meta['parallel']=dict(rank_ranges=ranges)
            actual=plane_statistics(paths,meta,1.4)[2]
            for key in expected:np.testing.assert_allclose(actual[key],expected[key],atol=1e-12)
            self.assertAlmostEqual(actual['F_uu'][-1],8/9)
            self.assertAlmostEqual(actual['R_uu'][-1],1)
            for files in (paths[:1],paths+paths[:1]):
                with self.assertRaises(ValueError):plane_statistics(files,meta,1.4)

    def test_direction_speed_and_cli(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);out=root/'analysis';meta=root/'meta.json'
            meta.write_text(json.dumps(dict(grid=[40,4,4],spacing=[1,1,1])))
            # Save two steps with physical times differing by 0.5; helper header step is patched below.
            for step,shock,time in [(10,20,1.),(20,18,1.5)]:
                path=root/f'field_{step:06d}.slf';write_slf(path,state(shock=shock,direction=-1),time=time)
                import struct
                with path.open('r+b') as f:f.seek(36);f.write(struct.pack('<i',step))
            args=[str(root),'--output-dir',str(out),'--meta',str(meta),'--gamma','1.4',
                  '--direction','-1','--search','10','30','--upstream','2','8','--downstream','2','8']
            self.assertEqual(main(args),0)
            import csv
            with (out/'shock_history.csv').open() as f:rows=list(csv.DictReader(f))
            self.assertAlmostEqual(float(rows[0]['shock_speed']),-4)
            self.assertEqual(json.loads((out/'analysis.json').read_text())['status'],'complete')
            with np.load(out/'shock_surface_00000020.npz') as saved:
                np.testing.assert_array_equal(saved['shock_x'],np.full((4,4),18.))
                self.assertEqual(float(saved['time']),1.5)
                self.assertTrue(saved['valid'].all())
                np.testing.assert_array_equal(saved['y'],np.arange(4)+.5)
            self.assertEqual(float(rows[1]['surface_x_std']),0)
            self.assertEqual(main(args),1)

    def test_tracking_zero_turbulence_and_invalid_data(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'field_000010.slf';q=state()
            q[...,4]-=.5*q[...,1]**2/q[...,0];q[...,1]=0
            write_slf(path,q)
            meta=dict(grid=[40,4,4],spacing=[1,1,1])
            prof=plane_statistics([path],meta,1.4)[2]
            result=diagnose(prof,1,[10,30],[2,8],[2,8])
            self.assertTrue(np.isnan(result['amplification_tke']))
            with self.assertRaises(ValueError):diagnose(prof,1,[10,30],[2,8],[2,8],previous=10,max_shift=2)
            q[1,1,1,0]=-1;write_slf(path,q)
            with self.assertRaises(ValueError):plane_statistics([path],meta,1.4)


if __name__=='__main__':unittest.main()
