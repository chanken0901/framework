"""Fortran executable comparisons. Set RF_FORTRAN_BUILD to the CMake build directory."""
import os
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'reference/python'))
sys.path.insert(0, str(ROOT/'tools'))
from export_mechanism import export
from reactingflow.importer import import_cantera


class FortranTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get('RF_FORTRAN_BUILD'):
            raise unittest.SkipTest('Set RF_FORTRAN_BUILD to enable compiled Fortran tests')
        import cantera as ct
        import numpy as np
        cls.ct, cls.np = ct, np
        suffix = '.exe' if os.name == 'nt' else ''
        cls.probe = Path(os.environ['RF_FORTRAN_BUILD'])/('rf_probe'+suffix)
        cls.reactor = cls.probe.with_name('rf_reactor'+suffix)
        if not cls.probe.is_file() or not cls.reactor.is_file():
            raise RuntimeError('Build rf_probe and rf_reactor first')
        cls.hydrogen = Path(ct.__file__).parent/'data/h2o2.yaml'

    def compare(self, ref, imported, path, t, p, y=None):
        y = [1/ref.n_species]*ref.n_species if y is None else y
        ref.TPY = t, p, y
        text = f'{t} {p}\n'+' '.join(map(str, y))+'\n'
        run = subprocess.run([str(self.probe), str(path)], input=text, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr+run.stdout)
        rows = [self.np.fromstring(line, sep=' ') for line in run.stdout.splitlines()]
        expect = [ref.cp_mass, ref.cv_mass, ref.enthalpy_mass, ref.int_energy_mass,
                  imported.gas.properties(t, y, p)['gas_constant'], ref.entropy_mass, ref.density, t,
                  ref.heat_release_rate]
        self.np.testing.assert_allclose(rows[0][:7], expect[:7], rtol=3.e-10, atol=1.e-6)
        self.assertAlmostEqual(rows[0][7], t, delta=1.e-5)  # NASA switch discontinuity
        self.np.testing.assert_allclose(rows[0][8], expect[8], rtol=3.e-10, atol=1.e-6)
        for a,b in zip(rows[1:], [ref.forward_rates_of_progress*1000, ref.reverse_rates_of_progress*1000,
                                 ref.net_rates_of_progress*1000, ref.net_production_rates*ref.molecular_weights]):
            self.np.testing.assert_allclose(a,b,rtol=3.e-9,atol=1.e-10)

    def test_nasa7_hydrogen_and_gri30(self):
        for name in ('h2o2.yaml','gri30.yaml'):
            path = self.hydrogen.with_name(name)
            own = import_cantera(path)
            ref = self.ct.Solution(str(path))
            with tempfile.TemporaryDirectory() as tmp:
                data = Path(tmp)/'mechanism.rf'
                export(own,data)
                for t in (400.,900.,1500.,2800.):
                    for p in (1.e3,1.e5,1.e7):
                        self.compare(ref,own,data,t,p)
                self.compare(ref,own,data,1000.,1.e5,[1.]+[0.]*(ref.n_species-1))

    def test_nasa9(self):
        ct = self.ct
        species = [s for s in ct.Species.list_from_file('airNASA9.yaml') if s.name in ('N2','O2','NO')]
        ref = ct.Solution(thermo='ideal-gas',kinetics='gas',species=species,reactions=[])
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'source.yaml'; data = Path(tmp)/'mechanism.rf'
            ref.write_yaml(source); own=import_cantera(source); export(own,data)
            for t in (350.,800.,1500.,2800.): self.compare(ref,own,data,t,230000.)

    def test_pressure_dependent_rates(self):
        ct = self.ct
        kinds = [(ct.LindemannRate,()),(ct.TroeRate,(.5,100,1000)),
                 (ct.TroeRate,(.5,100,1000,5000)),(ct.SriRate,(1.2,100,1500,.9,.1))]
        rates = [(cls(low=ct.ArrheniusRate(1e12,.1,2e6),high=ct.ArrheniusRate(1e8,-.2,1e6),
                      falloff_coeffs=params), True) for cls,params in kinds]
        rates += [(ct.PlogRate([(1e4,ct.ArrheniusRate(1e8,0,1e6)),
                               (1e4,ct.ArrheniusRate(-1e7,0,1e6)),
                               (1e6,ct.ArrheniusRate(2e9,.1,2e6))]),False),
                  (ct.ChebyshevRate(temperature_range=(300,3000),pressure_range=(1e3,1e7),
                                   data=[[8,.2,-.02],[.5,.1,.04],[.03,.01,.005]]),False)]
        for rate,falloff in rates:
            reaction=ct.Reaction(equation='H + O2 (+M) <=> HO2 (+M)' if falloff else 'H + O2 => HO2',rate=rate)
            if falloff: reaction.third_body.efficiencies={'H2O':6,'AR':.7}
            ref=ct.Solution(thermo='ideal-gas',kinetics='gas',species=ct.Solution(str(self.hydrogen)).species(),
                            reactions=[reaction])
            with tempfile.TemporaryDirectory() as tmp:
                source=Path(tmp)/'source.yaml'; data=Path(tmp)/'mechanism.rf'
                ref.write_yaml(source); own=import_cantera(source); export(own,data)
                pressures=(1e3,1e5,1e7) if isinstance(rate,ct.ChebyshevRate) else (1e2,1e5,1e8)
                for t in (500.,1000.,2500.):
                    for p in pressures: self.compare(ref,own,data,t,p)

    def run_reactor(self, mode, tol, initial_temperature=1000.):
        ct=self.ct
        ref=ct.Solution(str(self.hydrogen)); ref.TPX=initial_temperature,101325.,'H2:2,O2:1,N2:3.76'
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp); data=tmp/'mechanism.rf'; inp=tmp/'reactor.in'; out=tmp/'history.csv'
            export(import_cantera(self.hydrogen),data)
            inp.write_text(f"&reactor mode='{mode}', temperature={initial_temperature}, rtol={tol}, "
                           f"atol_species={tol*1e-7}, atol_temperature={tol*10}, end_time=0.001 /\n"+
                           ' '.join(map(str,ref.Y))+'\n',encoding='ascii')
            run=subprocess.run([str(self.reactor),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            text=out.read_text()
            self.assertTrue(text.rstrip().endswith('# SUCCESS'))
            lines=[line for line in text.splitlines() if not line.startswith('#')]
            values=self.np.loadtxt(lines[1:],delimiter=',',ndmin=2)
            diagnostics=dict(line[2:].split('=',1) for line in text.splitlines() if line.startswith('# ') and '=' in line)
            again=subprocess.run([str(self.reactor),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertNotEqual(again.returncode,0)
            self.assertEqual(out.read_text(),text)
        kind=ct.IdealGasReactor if mode=='constant_volume' else ct.IdealGasConstPressureReactor
        reactor=kind(ref,clone=True); network=ct.ReactorNet([reactor]); network.rtol=1.e-12;network.atol=1.e-18
        expected=[]
        for row in values:
            network.advance(row[0]);expected.append([reactor.T,reactor.phase.P,reactor.density,*reactor.phase.Y])
        return values,self.np.asarray(expected),diagnostics

    def test_reactor_histories_and_conservation(self):
        for mode in ('constant_volume','constant_pressure'):
            values,ref,diagnostics=self.run_reactor(mode,1e-10)
            self.np.testing.assert_allclose(values[:,1:4],ref[:,:3],rtol=3e-6)
            self.np.testing.assert_allclose(values[:,4:],ref[:,3:],rtol=5e-4,atol=3e-7)
            self.assertEqual(values[-1,0],.001)
            self.assertGreaterEqual(values[:,4:].min(),0)
            self.assertLess(float(diagnostics['energy_relative_error']),2e-7)
            self.assertLess(float(diagnostics['element_relative_error']),1e-12)
            self.assertLess(float(diagnostics['mass_sum_error']),1e-12)
            self.assertGreater(float(diagnostics['ignition_delay_s']),0)

    def test_tolerance_convergence(self):
        for mode in ('constant_volume','constant_pressure'):
            loose,lref,ld=self.run_reactor(mode,1e-5,1100.)
            tight,tref,td=self.run_reactor(mode,1e-9,1100.)
            self.assertLess(abs(tight[:,1]-tref[:,0]).max(),abs(loose[:,1]-lref[:,0]).max()*.1)
            self.assertLess(float(td['energy_relative_error']),float(ld['energy_relative_error'])*.1)

    def test_invalid_composition_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            data=Path(tmp)/'mechanism.rf'; export(import_cantera(self.hydrogen),data)
            for y in ([.2]*10,[-.1]+[1.1]+[0.]*8):
                run=subprocess.run([str(self.probe),str(data)],input='1000 101325\n'+' '.join(map(str,y)),
                                   capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)

    def test_custom_order_zero_collider_and_equilibrium(self):
        ct=self.ct
        species=ct.Solution(str(self.hydrogen)).species()
        custom=ct.Reaction(equation='H + O2 => HO2',rate=ct.ArrheniusRate(1e8,.1,2e6))
        custom.allow_nonreactant_orders=True
        custom.orders={'H':.5,'O2':1.2,'H2O':.1}
        zero=ct.Reaction(equation='H + O2 (+M) <=> HO2 (+M)',
                        rate=ct.LindemannRate(low=ct.ArrheniusRate(1e12,0,0),high=ct.ArrheniusRate(1e8,0,0)))
        zero.third_body.default_efficiency=0
        zero.third_body.efficiencies={'AR':1}
        for reaction in (custom,zero):
            ref=ct.Solution(thermo='ideal-gas',kinetics='gas',species=species,reactions=[reaction])
            with tempfile.TemporaryDirectory() as tmp:
                source=Path(tmp)/'source.yaml'; data=Path(tmp)/'mechanism.rf'
                ref.write_yaml(source); own=import_cantera(source); export(own,data)
                if reaction is zero:
                    ref.TPX=1000,1e5,'H:1,O2:1,HO2:1'
                    self.compare(ref,own,data,1000.,1e5,ref.Y.tolist())
                else: self.compare(ref,own,data,1000.,1e5)
        ref=ct.Solution(str(self.hydrogen)); ref.TPX=1500,1e5,'H2:2,O2:1,N2:4'; ref.equilibrate('TP')
        with tempfile.TemporaryDirectory() as tmp:
            data=Path(tmp)/'mechanism.rf'; export(import_cantera(self.hydrogen),data)
            run=subprocess.run([str(self.probe),str(data)],input='1500 100000\n'+' '.join(map(str,ref.Y)),
                               capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            rows=[self.np.fromstring(line,sep=' ') for line in run.stdout.splitlines()]
            self.assertLess(abs(rows[3]).max()/rows[1].max(),1e-10)


if __name__=='__main__': unittest.main()
