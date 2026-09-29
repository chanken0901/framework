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

    def test_frozen_shock_reference(self):
        from scipy.optimize import brentq
        gas = self.ct.Solution(str(self.hydrogen))
        exe = self.probe.with_name('rf_normal_shock'+self.probe.suffix)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            data, inp, out = tmp/'mechanism.rf', tmp/'shock.in', tmp/'shock.csv'
            export(import_cantera(self.hydrogen), data)
            for mach in (1.1, 2., 3.):
                gas.TPX = 300., 101325., 'H2:2,O2:1,N2:3.76'
                y = gas.Y.copy()
                rho0, h0, speed = gas.density, gas.enthalpy_mass, mach*gas.sound_speed
                def balance(ratio):
                    p = 101325.+rho0*speed**2*(1-1/ratio)
                    gas.DPY = rho0*ratio, p, y
                    return gas.enthalpy_mass+.5*(speed/ratio)**2-h0-.5*speed**2
                ratio = brentq(balance, 1.00001, rho0*speed**2/101325.)
                balance(ratio)
                expected = [speed, gas.T, gas.P, gas.density, speed*(1-1/ratio)]
                inp.write_text(f'&shock temperature=300, pressure=101325, mach={mach} /\n'
                               +' '.join(map(str,y))+'\n')
                run = subprocess.run([str(exe),str(data),str(inp),str(out)], capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                rows = [s for s in out.read_text().splitlines() if not s.startswith('#')]
                self.np.testing.assert_allclose(self.np.fromstring(rows[1],sep=','),expected,rtol=1.e-8)
                again = subprocess.run([str(exe),str(data),str(inp),str(out)], capture_output=True,text=True)
                self.assertNotEqual(again.returncode,0)
                out.unlink()
            for mach in (.9, 100.):
                inp.write_text(f'&shock mach={mach} /\n'+' '.join(map(str,y))+'\n')
                run = subprocess.run([str(exe),str(data),str(inp),str(out)], capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)
                self.assertFalse(out.exists())

    def test_znd_cantera_and_inert(self):
        from scipy.optimize import brentq
        from scipy.integrate import solve_ivp
        np = self.np
        gas = self.ct.Solution(str(self.hydrogen))
        exe = self.probe.with_name('rf_znd'+self.probe.suffix)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            data, inp, out = tmp/'mechanism.rf', tmp/'znd.in', tmp/'znd.csv'
            export(import_cantera(self.hydrogen), data)
            for composition, mach, tend in [('N2:1',3.,1.e-6),('H2:2,O2:1,N2:3.76',5.,1.e-6),
                                             ('H2:2,O2:1,N2:3.76',5.,1.e-4)]:
                gas.TPX = 300.,101325.,composition
                y0 = gas.Y.copy()
                rho0,h0,speed = gas.density,gas.enthalpy_mass,mach*gas.sound_speed
                flux = rho0*speed
                def balance(ratio):
                    gas.DPY = rho0*ratio,101325+rho0*speed**2*(1-1/ratio),y0
                    return gas.enthalpy_mass+.5*(speed/ratio)**2-h0-.5*speed**2
                ratio = brentq(balance,1.00001,rho0*speed**2/101325.)
                balance(ratio)
                initial = np.r_[gas.P,gas.density,0.,y0]
                def rhs(time,state):
                    p,rho,x = state[:3]
                    gas.DPY = rho,p,state[3:]
                    u = flux/rho
                    dy = gas.net_production_rates*gas.molecular_weights/rho
                    sigma = np.dot(gas.mean_molecular_weight/gas.molecular_weights
                                   -gas.partial_molar_enthalpies/gas.molecular_weights/(gas.cp_mass*gas.T),dy)
                    drho = -rho*sigma/(1-u*u/gas.sound_speed**2)
                    return np.r_[u*u*drho,drho,u,dy]
                ref = solve_ivp(rhs,[0,tend],initial,method='Radau',rtol=1.e-10,atol=1.e-15,
                                max_step=tend/100,dense_output=True)
                self.assertTrue(ref.success,ref.message)
                for tol in (1.e-8,1.e-10):
                    eqtol=1.e-6 if tend>1.e-5 else 0.
                    maxstep=tend/(100 if tol==1.e-8 else 200)
                    inp.write_text(f'&znd mach={mach}, end_time={tend}, rtol={tol}, max_step={maxstep}, '
                                   f'equilibrium_tolerance={eqtol} /\n'
                                   +' '.join(map(str,y0))+'\n')
                    run = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                    self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                    self.assertIn('# SUCCESS',out.read_text())
                    rows = np.loadtxt([line for line in out.read_text().splitlines()
                                       if not line.startswith(('#','residence_time'))],delimiter=',')
                    for index in sorted(set(range(0,len(rows),max(1,len(rows)//10))) | {len(rows)-1}):
                        row = rows[index]
                        expected = ref.sol(row[0])
                        gas.DPY = expected[1],expected[0],expected[3:]
                        np.testing.assert_allclose(row[1:6],
                            [expected[2],gas.T,expected[0],expected[1],flux/expected[1]],rtol=2.e-5,atol=1.e-10)
                        np.testing.assert_allclose(row[7:],expected[3:],rtol=2.e-4,atol=2.e-8)
                    if composition=='N2:1':
                        np.testing.assert_allclose(rows[-1,2:6],rows[0,2:6],rtol=1.e-10)
                    else:
                        self.assertGreater(rows[-1,2]-rows[0,2],100.)
                    np.testing.assert_allclose(rows[:,4]*rows[:,5],flux,rtol=1.e-12)
                    np.testing.assert_allclose(rows[:,3]+rows[:,4]*rows[:,5]**2,
                                               101325+rho0*speed**2,rtol=1.e-12)
                    saved = out.read_bytes()
                    again = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                    self.assertNotEqual(again.returncode,0)
                    self.assertEqual(saved,out.read_bytes())
                    out.unlink()
                inp.write_text('&znd sonic_margin=0.99 /\n'+' '.join(map(str,y0))+'\n')
                run = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)
                self.assertNotIn('# SUCCESS',out.read_text())
                out.unlink()
            inp.write_text("&znd speed_mode='cj', overdrive=1.05, end_time=1.e-4, max_step=1.e-7, "
                           "equilibrium_tolerance=1.e-6 /\n"+' '.join(map(str,y0))+'\n')
            run = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertIn('# equilibrium_max_abs_Y_error=',out.read_text())
            self.assertIn('# SUCCESS',out.read_text())
            out.unlink()
            example = ROOT/'examples/znd_h2_air.in'
            run = subprocess.run([str(exe),str(data),str(example),str(out)],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            self.assertIn('# SUCCESS',out.read_text())
            out.unlink()
            inp.write_text('&znd max_steps=1 /\n'+' '.join(map(str,y0))+'\n')
            run = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            self.assertNotIn('# SUCCESS',out.read_text())
            out.unlink()
            inp.write_text('&znd mach=5, end_time=1.e-8, equilibrium_tolerance=1.e-6 /\n'
                           +' '.join(map(str,y0))+'\n')
            run=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            self.assertNotIn('# SUCCESS',out.read_text())

    def test_cj_equilibrium_reference(self):
        from scipy.optimize import brentq, minimize_scalar
        gas = self.ct.Solution(str(self.hydrogen))
        exe = self.probe.with_name('rf_cj'+self.probe.suffix)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            data, inp, out = tmp/'mechanism.rf',tmp/'cj.in',tmp/'cj.csv'
            export(import_cantera(self.hydrogen),data)
            for p0 in (101325.,202650.):
                gas.TPX = 300.,p0,'H2:2,O2:1,N2:3.76'
                y0,rho0,e0 = gas.Y.copy(),gas.density,gas.int_energy_mass
                def wave(ratio):
                    def energy(t):
                        gas.TDY = t,rho0*ratio,y0
                        gas.equilibrate('TV',solver='gibbs',rtol=1.e-11,max_steps=2000)
                        return gas.int_energy_mass-e0-.5*(gas.P+p0)*(1/rho0-1/(rho0*ratio))
                    temp = brentq(energy,500.,3500.,xtol=1.e-7)
                    energy(temp)
                    return (gas.P-p0)/(rho0*(1-1/ratio))
                ref = minimize_scalar(wave,bounds=(1.2,2.5),method='bounded',options={'xatol':1.e-8})
                self.assertTrue(ref.success)
                d2 = wave(ref.x)
                expected = self.np.r_[d2**.5,gas.T,gas.P,gas.density,gas.Y]
                inp.write_text(f'&cj temperature=300, pressure={p0} /\n'+' '.join(map(str,y0))+'\n')
                run = subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                row = self.np.loadtxt([s for s in out.read_text().splitlines()
                                      if not s.startswith(('#','cj_speed'))],delimiter=',')
                self.np.testing.assert_allclose(row[:4],expected[:4],rtol=2.e-6)
                self.np.testing.assert_allclose(row[4:],expected[4:],rtol=1.e-4,atol=1.e-9)
                # Independent equilibrium isentropic sound speed (not frozen sound speed).
                gas.TDY = row[1],row[3],row[4:]
                entropy=gas.entropy_mass
                pressures=[]
                for factor in (1-1.e-4,1+1.e-4):
                    gas.TDY=row[1],row[3],row[4:]
                    gas.SV=entropy,1/(row[3]*factor)
                    gas.equilibrate('SV',rtol=1.e-10)
                    pressures.append(gas.P)
                aeq=((pressures[1]-pressures[0])/(2.e-4*row[3]))**.5
                self.assertAlmostEqual(row[0]*rho0/row[3]/aeq,1.,delta=2.e-5)
                saved=out.read_bytes()
                again=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
                self.assertNotEqual(again.returncode,0)
                self.assertEqual(saved,out.read_bytes())
                out.unlink()
            inp.write_text('&cj ratio_min=1.01, ratio_max=1.1 /\n'+' '.join(map(str,y0))+'\n')
            run=subprocess.run([str(exe),str(data),str(inp),str(out)],capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            self.assertFalse(out.exists())

    def test_equilibrium_tv(self):
        exe=self.probe.with_name('rf_equilibrium_probe'+self.probe.suffix)
        for mechanism,composition in [('h2o2.yaml','H2:2,O2:1,N2:3.76'),
                                      ('h2o2.yaml','H2:2,O2:1'),
                                      ('gri30.yaml','CH4:1,O2:2,N2:7.52')]:
            source=self.hydrogen.with_name(mechanism)
            gas=self.ct.Solution(str(source))
            with tempfile.TemporaryDirectory() as tmp:
                data=Path(tmp)/'mechanism.rf';export(import_cantera(source),data)
                for temp in (500.,1500.,3000.):
                    for rho in (.1,3.):
                        gas.TPX=300.,101325.,composition
                        y0=gas.Y.copy()
                        gas.TDY=temp,rho,y0;gas.equilibrate('TV',rtol=1.e-11,max_steps=3000)
                        run=subprocess.run([str(exe),str(data)],input=f'{temp} {rho}\n'+' '.join(map(str,y0))+'\n',
                                           capture_output=True,text=True)
                        self.assertEqual(run.returncode,0,run.stdout+run.stderr)
                        actual=self.np.fromstring(run.stdout,sep=' ')
                        self.np.testing.assert_allclose(actual,gas.Y,rtol=2.e-5,atol=1.e-9)

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
