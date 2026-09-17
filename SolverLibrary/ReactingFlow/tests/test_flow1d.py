"""R4a executable validation; Python is only a test harness/reference, not the solver."""
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference/python'))
sys.path.insert(0,str(ROOT/'tools'))
from export_mechanism import export
from reactingflow.importer import import_cantera


class FlowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get('RF_FORTRAN_BUILD'):
            raise unittest.SkipTest('Set RF_FORTRAN_BUILD for compiled CFD tests')
        import cantera as ct
        import numpy as np
        cls.ct,cls.np=ct,np
        cls.exe=Path(os.environ['RF_FORTRAN_BUILD'])/('rf_flow1d.exe' if os.name=='nt' else 'rf_flow1d')
        if not cls.exe.is_file(): raise RuntimeError('Build rf_flow1d first')
        cls.h2=Path(ct.__file__).parent/'data/h2o2.yaml'

    def run_flow(self,source,controls,yl,yr,success=True):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);data=tmp/'mechanism.rf';inp=tmp/'flow.in';out=tmp/'flow.csv'
            export(import_cantera(source),data)
            inp.write_text('&flow1d '+controls+' /\n'+' '.join(map(str,yl))+'\n'+' '.join(map(str,yr))+'\n')
            run=subprocess.run([str(self.exe),str(data),str(inp),str(out)],capture_output=True,text=True,timeout=120)
            if not success:
                self.assertNotEqual(run.returncode,0)
                if out.exists(): self.assertNotIn('# SUCCESS',out.read_text())
                return
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            text=out.read_text()
            self.assertTrue(text.rstrip().endswith('# SUCCESS'))
            rows=[line for line in text.splitlines() if not line.startswith('#')]
            values=self.np.loadtxt(rows[1:],delimiter=',',ndmin=2)
            diagnostics=dict(line[2:].split('=',1) for line in text.splitlines() if line.startswith('# ') and '=' in line)
            self.assertGreaterEqual(int(diagnostics['rejected_steps']),0)
            return values,diagnostics

    def test_dirichlet_uniform_and_reactive_conservation(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76'
        y=gas.Y.tolist();fractions=','.join(map(str,y))
        for chemistry in (False,True):
            with self.subTest(chemistry=chemistry):
                controls=("nx=8,length=0.1,interface_x=0.05,end_time=0.000002,max_dt=0.0000001,"
                    "left_temperature=1100,right_temperature=1100,left_velocity=20,right_velocity=20,"
                    "left_bc='dirichlet',right_bc='dirichlet',reconstruction='muscl',"
                    "transport_model='constant',viscosity=0.00002,thermal_conductivity=0.03,mass_diffusivity=0.00001,"
                    f"chemistry={'.true.' if chemistry else '.false.'},"
                    f"left_boundary_temperature=1100,left_boundary_pressure=101325,left_boundary_velocity=20,left_boundary_y={fractions},"
                    f"right_boundary_temperature=1100,right_boundary_pressure=101325,right_boundary_velocity=20,right_boundary_y={fractions}")
                values,d=self.run_flow(self.h2,controls,y,y)
                for key in ('mass_error','momentum_error','energy_error','element_error'):
                    self.assertLess(float(d[key]),1e-8)
                self.assertIn('fixed_state_left',d);self.assertIn('fixed_state_right',d)
                final=values[values[:,0]==values[-1,0]]
                self.assertGreater(final[:,6].min(),0)
                if not chemistry:
                    self.np.testing.assert_allclose(final[:,3:],values[values[:,0]==0][:,3:],rtol=1e-10,atol=1e-12)

    def test_dirichlet_missing_or_unused_state_rejected(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        for controls in ("left_bc='dirichlet'", "left_boundary_temperature=1100",
                         "left_bc='dirichlet',left_boundary_temperature=1100,left_boundary_pressure=101325"):
            with self.subTest(controls=controls):
                self.run_flow(self.h2,controls,y,y,success=False)

    def test_dirichlet_example_and_right_inflow(self):
        example=(ROOT/'examples/flow1d_h2_inflow.in').read_text()
        body=example.split('&flow1d',1)[1]
        controls,rows=body.split('/',1)
        yl,yr=([float(v) for v in line.split()] for line in rows.strip().splitlines())
        for right in (False,True):
            selected=controls
            if right:
                selected=selected.replace("left_bc='dirichlet', right_bc='outflow'",
                    "left_bc='outflow', right_bc='dirichlet'").replace('left_boundary_', 'right_boundary_')
                selected=selected.replace('velocity=20', 'velocity=-20')
            values,d=self.run_flow(self.h2,selected,yl,yr)
            self.assertEqual(values[-1,1],0.000002)
            for key in ('mass_error','momentum_error','energy_error','element_error'):
                self.assertLess(float(d[key]),1e-8)

    def test_uniform_reactive_matches_cantera(self):
        ct=self.ct
        gas=ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76'
        y=gas.Y.tolist()
        values,d=self.run_flow(self.h2,"nx=4,length=100,interface_x=50,end_time=0.0002,max_dt=0.00001,"+
            "write_every=1,chemistry=.true.,left_bc='periodic',right_bc='periodic',"+
            "left_temperature=1100,right_temperature=1100",y,y)
        reactor=ct.IdealGasReactor(gas,clone=True);net=ct.ReactorNet([reactor]);net.rtol=1e-12;net.atol=1e-18
        for time in self.np.unique(values[:,1]):
            net.advance(time)
            row=values[values[:,1]==time]
            self.np.testing.assert_allclose(row[:,5],reactor.T,rtol=3e-6)
            self.np.testing.assert_allclose(row[:,6],reactor.phase.P,rtol=3e-6)
            self.np.testing.assert_allclose(row[:,7:],self.np.tile(reactor.phase.Y,(4,1)),rtol=5e-4,atol=1e-7)
        self.assertLess(float(d['energy_error']),1e-12)
        self.assertEqual(int(d['rejected_steps']),0)
        self.assertLess(float(d['element_error']),1e-8)
        self.assertEqual(values[-1,1],.0002)

    def test_periodic_contact_conserves_species(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';yl=gas.Y.tolist()
        gas.TPX=1100,101325,'AR:1';yr=gas.Y.tolist()
        values,d=self.run_flow(self.h2,"nx=24,length=1,end_time=0.0001,max_dt=0.00001,"+
            "left_bc='periodic',right_bc='periodic',left_velocity=100,right_velocity=100",yl,yr)
        first=values[values[:,0]==0];last=values[values[:,0]==values[-1,0]]
        self.np.testing.assert_allclose((first[:,3,None]*first[:,7:]).sum(0),
                                        (last[:,3,None]*last[:,7:]).sum(0),rtol=1e-12,atol=1e-13)
        for key in ('mass_error','momentum_error','energy_error','element_error'):
            self.assertLess(float(d[key]),1e-10)

    def test_sod_grid_convergence(self):
        ct=self.ct
        species=ct.Species('AR',{'Ar':1})
        species.thermo=ct.NasaPoly2(1e-4,1e5,101325,[1000,3.5,0,0,0,0,0,0,3.5,0,0,0,0,0,0])
        gas=ct.Solution(thermo='ideal-gas',kinetics='gas',species=[species],reactions=[])
        rmix=ct.gas_constant/gas.mean_molecular_weight
        errors=[]
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'ideal.yaml';gas.write_yaml(source)
            for nx in (80,160):
                values,d=self.run_flow(source,f'nx={nx},end_time=0.2,max_dt=0.01,write_every=100000,'+
                    f'left_temperature={1/rmix},right_temperature={.1/.125/rmix},left_pressure=1,right_pressure=0.1',
                    [1.],[1.])
                final=values[values[:,0]==values[-1,0]]
                exact=self.sod_density((final[:,2]-.5)/.2)
                errors.append(abs(final[:,3]-exact).mean())
                self.assertGreater(final[:,6].min(),0)
                self.assertLess(float(d['energy_error']),1e-10)
        self.assertLess(errors[1],errors[0]*.8)
        self.assertLess(errors[1],.035)

    def sod_density(self,xi):
        # Exact perfect-gas Sod solution, gamma=1.4; pressure root independent of FV code.
        import math
        gamma=1.4;al=math.sqrt(gamma);ar=math.sqrt(gamma*.1/.125)
        def fl(p): return 2*al/(gamma-1)*(p**((gamma-1)/(2*gamma))-1)
        def fr(p): return (p-.1)*math.sqrt(2/((gamma+1)*.125)/(p+(gamma-1)/(gamma+1)*.1))
        lo,hi=.1,1.
        for _ in range(80):
            p=(lo+hi)/2
            if fl(p)+fr(p)>0: hi=p
            else: lo=p
        p=(lo+hi)/2;u=(fr(p)-fl(p))/2
        astar=al*p**((gamma-1)/(2*gamma));tail=u-astar
        shock=ar*math.sqrt((gamma+1)/(2*gamma)*p/.1+(gamma-1)/(2*gamma))
        alpha=(gamma-1)/(gamma+1);right=.125*(p/.1+alpha)/(alpha*p/.1+1)
        out=[]
        for v in xi:
            if v < -al: rho=1.
            elif v < tail: rho=(2/(gamma+1)*(al-(gamma-1)*v/2)/al)**(2/(gamma-1))
            elif v < u: rho=p**(1/gamma)
            elif v < shock: rho=right
            else: rho=.125
            out.append(rho)
        return self.np.asarray(out)

    def test_bad_configuration_and_limit_rejected(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        for controls in ["cfl=2", "left_bc='periodic'", "max_steps=1,end_time=1", "chemistry_rtol=0,chemistry=.true."]:
            self.run_flow(self.h2,controls,y,y,success=False)

    def test_nonuniform_reacting_time_refinement(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76';y=gas.Y.tolist()
        solutions=[]
        for dt in (1.e-5,5.e-6,2.5e-6):
            values,d=self.run_flow(self.h2,f'nx=8,length=1,end_time=0.00004,max_dt={dt},write_every=100000,'+
                "chemistry=.true.,left_bc='reflecting',right_bc='reflecting',"+
                "left_temperature=1300,right_temperature=1100",y,y)
            solutions.append(values[values[:,0]==values[-1,0]])
            self.assertLess(float(d['mass_error']),1e-12)
            self.assertLess(float(d['energy_error']),1e-12)
            self.assertLess(float(d['element_error']),1e-8)
        # Temporal refinement on the same grid; not a spatial accuracy/detonation claim.
        coarse=abs(solutions[0][:,5]-solutions[2][:,5]).max()
        medium=abs(solutions[1][:,5]-solutions[2][:,5]).max()
        self.assertGreater(coarse,1.e-6)
        self.assertLess(medium,coarse*.6)

    def test_zero_transport_preserves_euler_result(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        controls="nx=12,end_time=0.00004,left_temperature=1300,right_temperature=1100,write_every=1"
        baseline,_=self.run_flow(self.h2,controls,y,y)
        explicit,_=self.run_flow(self.h2,controls+",transport_model='constant'",y,y)
        self.np.testing.assert_array_equal(explicit,baseline)

    def test_reacting_transport_conservation_and_effect(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76';y=gas.Y.tolist()
        controls="nx=8,length=1,end_time=0.00004,max_dt=0.0000025,write_every=100000,"+\
            "chemistry=.true.,left_bc='reflecting',right_bc='reflecting',"+\
            "left_temperature=1300,right_temperature=1100"
        base,_=self.run_flow(self.h2,controls,y,y)
        values,d=self.run_flow(self.h2,controls+",transport_model='constant',viscosity=0.1,"+\
            "bulk_viscosity=0.05,thermal_conductivity=100,mass_diffusivity=0.1",y,y)
        final=values[values[:,0]==values[-1,0]]
        ref=base[base[:,0]==base[-1,0]]
        self.assertGreater(abs(final[:,5]-ref[:,5]).max(),1e-5)
        self.assertGreater(final[:,3].min(),0)
        self.assertGreater(final[:,6].min(),0)
        self.assertGreaterEqual(final[:,7:].min(),0)
        for key in ('mass_error','momentum_error','energy_error','element_error'):
            self.assertLess(float(d[key]),1e-8)
        self.assertEqual(d['transport_model'],'constant')
        self.assertEqual(values[-1,1],.00004)

    def test_reconstruction_configuration_and_default(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        controls="nx=12,end_time=0.00004,left_temperature=1300,right_temperature=1100,write_every=1"
        baseline,d=self.run_flow(self.h2,controls,y,y)
        explicit,_=self.run_flow(self.h2,controls+",reconstruction='first_order'",y,y)
        self.np.testing.assert_array_equal(explicit,baseline)
        self.assertEqual(d['reconstruction'],'first_order')
        self.run_flow(self.h2,"reconstruction='unknown'",y,y,success=False)

    def test_muscl_reacting_transport_conservation(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76';y=gas.Y.tolist()
        values,d=self.run_flow(self.h2,"nx=12,length=1,end_time=0.00004,max_dt=0.0000025,"+
            "write_every=100000,chemistry=.true.,reconstruction='muscl',"+
            "left_bc='reflecting',right_bc='reflecting',left_temperature=1300,right_temperature=1100,"+
            "transport_model='constant',viscosity=0.1,thermal_conductivity=100,mass_diffusivity=0.1",y,y)
        self.assertEqual(d['reconstruction'],'muscl')
        self.assertGreater(values[:,6].min(),0)
        self.assertGreaterEqual(values[:,7:].min(),0)
        self.np.testing.assert_allclose(values[:,7:].sum(1),1,atol=1e-12)
        for key in ('mass_error','momentum_error','energy_error','element_error'):
            self.assertLess(float(d[key]),1e-8)
        self.assertEqual(values[-1,1],.00004)

    def test_muscl_sod_accuracy(self):
        ct=self.ct
        species=ct.Species('AR',{'Ar':1})
        species.thermo=ct.NasaPoly2(1e-4,1e5,101325,[1000,3.5,0,0,0,0,0,0,3.5,0,0,0,0,0,0])
        gas=ct.Solution(thermo='ideal-gas',kinetics='gas',species=[species],reactions=[])
        rmix=ct.gas_constant/gas.mean_molecular_weight
        errors={}
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'ideal.yaml';gas.write_yaml(source)
            for method in ('first_order','muscl'):
                values,d=self.run_flow(source,f"nx=80,reconstruction='{method}',end_time=0.2,"+
                    "max_dt=0.01,write_every=100000,"+
                    f'left_temperature={1/rmix},right_temperature={.1/.125/rmix},left_pressure=1,right_pressure=0.1',
                    [1.],[1.])
                final=values[values[:,0]==values[-1,0]]
                errors[method]=abs(final[:,3]-self.sod_density((final[:,2]-.5)/.2)).mean()
                self.assertGreater(final[:,6].min(),0)
                self.assertLess(float(d['energy_error']),1e-10)
        self.assertLess(errors['muscl'],.7*errors['first_order'])

    def test_species_diffusion_equal_coefficients(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';yl=gas.Y.tolist()
        gas.TPX=1100,101325,'AR:1';yr=gas.Y.tolist()
        controls="nx=12,end_time=0.00004,write_every=1,reconstruction='muscl',"+\
            "left_bc='periodic',right_bc='periodic'"
        common,_=self.run_flow(self.h2,controls+",transport_model='constant',mass_diffusivity=0.1",yl,yr)
        individual,d=self.run_flow(self.h2,controls+",transport_model='species_constant',species_diffusivities=10*0.1",yl,yr)
        self.np.testing.assert_array_equal(individual,common)
        self.assertEqual(d['transport_model'],'species_constant')
        self.assertEqual(float(d['diffusivity_H2']),.1)

    def test_species_diffusion_reacting_conservation(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'H2:2,O2:1,N2:3.76';yl=gas.Y.tolist()
        gas.TPX=1100,101325,'H2:1,O2:1,N2:3.76';yr=gas.Y.tolist()
        controls="nx=12,length=1,end_time=0.00002,max_dt=0.000002,write_every=100000,"+\
            "reconstruction='muscl',chemistry=.true.,left_temperature=1300,right_temperature=1100,"+\
            "left_bc='reflecting',right_bc='reflecting',viscosity=0.1,thermal_conductivity=100"
        base,_=self.run_flow(self.h2,controls+",transport_model='constant',mass_diffusivity=0.1",yl,yr)
        values,d=self.run_flow(self.h2,controls+",transport_model='species_constant',"+\
            "species_diffusivities=0.3,0.4,0.2,0.1,0.2,0.15,0.1,0.1,0.08,0.08",yl,yr)
        final=values[values[:,0]==values[-1,0]]
        ref=base[base[:,0]==base[-1,0]]
        self.assertGreater(abs(final[:,7:]-ref[:,7:]).max(),1e-8)
        self.assertGreaterEqual(values[:,7:].min(),0)
        self.np.testing.assert_allclose(values[:,7:].sum(1),1,atol=1e-12)
        for key in ('mass_error','momentum_error','energy_error','element_error'):
            self.assertLess(float(d[key]),1e-8)

    def test_invalid_species_diffusion_rejected(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        for controls in ["transport_model='species_constant'",
                         "transport_model='species_constant',species_diffusivities=0.1,0.2",
                         "transport_model='species_constant',species_diffusivities=10*-0.1",
                         "transport_model='species_constant',species_diffusivities=10*0.1,mass_diffusivity=0.1",
                         "transport_model='constant',species_diffusivities=10*0.1",
                         "transport_model='none',species_diffusivities=10*0.1"]:
            self.run_flow(self.h2,controls,y,y,success=False)

    def test_invalid_transport_rejected(self):
        gas=self.ct.Solution(str(self.h2));gas.TPX=1100,101325,'N2:1';y=gas.Y.tolist()
        for controls in ["transport_model='unknown'", "viscosity=1", "mass_diffusivity=1",\
                         "transport_model='constant',viscosity=-1",\
                         "transport_model='constant',bulk_viscosity=-1",\
                         "transport_model='constant',thermal_conductivity=-1",\
                         "transport_model='constant',mass_diffusivity=-1"]:
            self.run_flow(self.h2,controls,y,y,success=False)


if __name__=='__main__': unittest.main()
