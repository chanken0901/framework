"""Independent Cantera transport comparisons and table-backed CFD regression."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from export_mechanism import export,import_cantera
from export_transport import export_transport
import cantera as ct
import numpy as np


class TableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get('RF_FORTRAN_BUILD'):
            raise unittest.SkipTest('Set RF_FORTRAN_BUILD')
        cls.build=Path(os.environ['RF_FORTRAN_BUILD'])
        cls.suffix='.exe' if os.name=='nt' else ''
        cls.source=Path(ct.__file__).parent/'data/h2o2.yaml'

    def fixture(self,tmp):
        data=tmp/'mechanism.rf';export(import_cantera(self.source),data)
        directory=tmp/'transport'
        export_transport(self.source,directory,np.geomspace(300,3500,201))
        return data,directory

    def probe(self,data,directory,t,p,y,grad,ok=True):
        text=f'{t} {p}\n'+' '.join(map(str,y))+'\n'+' '.join(map(str,grad))+'\n'
        run=subprocess.run([str(self.build/('rf_transport_probe'+self.suffix)),str(data),
                            str(directory/'properties.rf'),str(directory/'binary.rf')],input=text,capture_output=True,text=True)
        if not ok:
            self.assertNotEqual(run.returncode,0);return
        self.assertEqual(run.returncode,0,run.stdout+run.stderr)
        return [np.fromstring(line,sep=' ') for line in run.stdout.splitlines()]

    def test_properties_diffusion_and_flux(self):
        with tempfile.TemporaryDirectory() as td:
            data,directory=self.fixture(Path(td))
            gas=ct.Solution(str(self.source));gas.transport_model='mixture-averaged'
            for t in (300.,733.,1500.,3500.):
                for p in (50000.,202650.):
                    for x in (np.ones(gas.n_species),np.arange(1,gas.n_species+1),[1.]+[0.]*(gas.n_species-1)):
                        gas.TPX=t,p,x;y=gas.Y.copy()
                        grad=np.zeros(gas.n_species)
                        if np.all(y>0): grad[0]=1.e-4;grad[1]=-1.e-4
                        rows=self.probe(data,directory,t,p,y,grad)
                        np.testing.assert_allclose(rows[0][:2],[gas.viscosity,gas.thermal_conductivity],rtol=8.e-5)
                        self.assertTrue(np.all(rows[0][2:]>=rows[0][:2]))
                        expected=gas.binary_diff_coeffs.copy();np.fill_diagonal(expected,0)
                        np.testing.assert_allclose(rows[1:gas.n_species+1],expected,rtol=8.e-5,atol=1.e-30)
                        mw=gas.molecular_weights;mean=gas.mean_molecular_weight
                        gradx=mean/mw*grad-gas.X*mean*np.sum(grad/mw)
                        flux=-gas.density*mw/mean*gas.mix_diff_coeffs*gradx
                        flux-=y*np.sum(flux)
                        np.testing.assert_allclose(rows[-1][:-2],flux,rtol=1.e-4,atol=1.e-13)
                        self.assertAlmostEqual(sum(rows[-1][:-2]),0.,delta=1.e-18)
                        enthalpy=np.dot(gas.partial_molar_enthalpies/mw,flux)
                        np.testing.assert_allclose(rows[-1][-1],enthalpy,rtol=1.e-4,atol=1.e-8)

    def test_invalid_tables_and_export(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);data,directory=self.fixture(tmp)
            gas=ct.Solution(str(self.source));gas.TPX=1000,101325,'N2:1';y=gas.Y
            for t in (299.,3501.): self.probe(data,directory,t,101325,y,np.zeros(len(y)),False)
            path=directory/'properties.rf';original=path.read_text()
            for damaged in (original.replace('RF_TRANSPORT_TABLE_V1','BAD'),
                            original.replace('H2 ','UNKNOWN ',1),original+'EXTRA\n'):
                path.write_text(damaged)
                self.probe(data,directory,1000,101325,y,np.zeros(len(y)),False)
            path.write_text(original)
            with self.assertRaises(FileExistsError): export_transport(self.source,directory,[300.,3500.])
            for temps in ([100.,400.],[300.,300.],[300.,float('nan')]):
                with self.assertRaises(ValueError): export_transport(self.source,tmp/'bad',temps)
                self.assertFalse((tmp/'bad').exists())

    def test_temperature_grid_refinement(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);data=tmp/'mechanism.rf';export(import_cantera(self.source),data)
            gas=ct.Solution(str(self.source));gas.TPX=733.,101325.,np.ones(gas.n_species)
            reference=np.r_[gas.viscosity,gas.thermal_conductivity]
            errors=[]
            for count in (21,201):
                directory=tmp/str(count)
                export_transport(self.source,directory,np.geomspace(300,3500,count))
                rows=self.probe(data,directory,733.,101325.,gas.Y,np.zeros(gas.n_species))
                errors.append(np.max(np.abs(rows[0][:2]/reference-1)))
            self.assertLess(errors[1],errors[0]/10)

    def test_table_backed_reacting_flow(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=Path(td);data,directory=self.fixture(tmp)
            gas=ct.Solution(str(self.source));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76'
            row=' '.join(map(str,gas.Y))
            controls="nx=8,length=1,end_time=1.e-7,max_dt=1.e-8,chemistry=.true.,reconstruction='muscl'," \
                     "left_bc='periodic',right_bc='periodic',left_temperature=1100,right_temperature=1120," \
                     "left_velocity=1,right_velocity=-1,transport_model='mixture_averaged'," \
                     "viscosity_model='tabulated_wilke',conductivity_model='tabulated_mix'," \
                     "transport_file='transport/properties.rf',binary_diffusion_model='tabulated'," \
                     "binary_diffusion_file='transport/binary.rf'"
            inp=tmp/'flow.in';out=tmp/'flow.csv'
            inp.write_text('&flow1d '+controls+' /\n'+row+'\n'+row+'\n')
            exe=self.build/('rf_flow1d'+self.suffix)
            run=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertTrue(out.read_text().rstrip().endswith('# SUCCESS'))
            diag=dict(s[2:].split('=',1) for s in out.read_text().splitlines() if s.startswith('# ') and '=' in s)
            for key,value in diag.items():
                if key.endswith('_error'): self.assertLess(float(value),1.e-8)
            out.unlink()
            inp.write_text((ROOT/'examples/flow_tabulated_transport.in').read_text())
            run=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertTrue(out.read_text().rstrip().endswith('# SUCCESS'))
            for wrong in (controls.replace("conductivity_model='tabulated_mix'","conductivity_model='constant'"),
                          controls+',viscosity=1',controls.replace("transport_model='mixture_averaged'","transport_model='none'")):
                out.unlink();inp.write_text('&flow1d '+wrong+' /\n'+row+'\n'+row+'\n')
                run=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)
                # Keep a placeholder for the next iteration's cleanup.
                if not out.exists(): out.touch()
