from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import subprocess
import json
import math
import tempfile
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'reference/python'))
from reactingflow.importer import import_cantera
from reactingflow.mechanism import MechanismError
from reactingflow.reactor import integrate
from reactingflow.thermo import NASA, IdealGas
from reactingflow.kinetics import Kinetics, Reaction, Rate, Arrhenius


class ReactorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import scipy
            import numpy as np
            import cantera as ct
        except ImportError:
            raise unittest.SkipTest('SciPy and Cantera required for R3 reference tests')
        cls.ct, cls.np = ct, np
        cls.path = Path(ct.__file__).parent/'data/h2o2.yaml'
        cls.mechanism = import_cantera(cls.path)
        gas = ct.Solution(str(cls.path))
        gas.TPX = 1000., 101325., 'H2:2,O2:1,N2:3.76'
        cls.y = gas.Y.tolist()

    def run_own(self, mode='constant_volume', **kwargs):
        args = dict(temperature=1000., pressure=101325., mass_fractions=self.y,
                    end_time=.001, mode=mode)
        args.update(kwargs)
        return integrate(self.mechanism, **args)

    def reference(self, mode, times):
        ct = self.ct
        gas = ct.Solution(str(self.path))
        gas.TPY = 1000., 101325., self.y
        kind = ct.IdealGasReactor if mode == 'constant_volume' else ct.IdealGasConstPressureReactor
        reactor = kind(gas, clone=True)
        network = ct.ReactorNet([reactor])
        network.rtol, network.atol = 1.e-12, 1.e-18
        states = []
        for t in times:
            network.advance(t)
            states.append([reactor.T, reactor.phase.P, reactor.density, *reactor.phase.Y])
        return self.np.asarray(states)

    def test_ignition_trajectories_and_conservation(self):
        for mode in ('constant_volume', 'constant_pressure'):
            with self.subTest(mode=mode):
                result = self.run_own(mode, rtol=1.e-9, atol_species=1.e-16, atol_temperature=1.e-8)
                ref = self.reference(mode, result.time)
                self.np.testing.assert_allclose(result.temperature, ref[:, 0], rtol=2.e-6)
                self.np.testing.assert_allclose(result.pressure, ref[:, 1], rtol=3.e-6)
                self.np.testing.assert_allclose(result.density, ref[:, 2], rtol=3.e-6)
                self.np.testing.assert_allclose(result.mass_fractions, ref[:, 3:], rtol=5.e-4, atol=3.e-7)
                self.assertLess(result.diagnostics['mass_sum_error'], 1.e-12)
                self.assertLess(result.diagnostics['element_relative_error'], 1.e-12)
                # h2o2.yaml has a small NASA enthalpy discontinuity at initial T=1000 K.
                self.assertLess(result.diagnostics['energy_relative_error'], 2.e-7)
                self.assertTrue(all(min(y) >= 0 for y in result.mass_fractions))
                self.assertEqual(result.time[-1], .001)
                ignition_ref = self.reference(mode, [result.ignition_delay])[0, 0]
                self.assertAlmostEqual(ignition_ref, 1400., delta=.02)

    def test_tolerance_convergence(self):
        for mode in ('constant_volume', 'constant_pressure'):
            loose = self.run_own(mode, rtol=1.e-5, atol_species=1.e-12, atol_temperature=1.e-4)
            tight = self.run_own(mode, rtol=1.e-9, atol_species=1.e-16, atol_temperature=1.e-8)
            error = lambda r: abs(self.reference(mode, [r.ignition_delay])[0, 0]-1400.)
            self.assertLess(error(tight), error(loose)*.1)
            self.assertLess(tight.diagnostics['energy_relative_error'], loose.diagnostics['energy_relative_error'])

    def test_nasa_switch_energy_floor_is_identified(self):
        gas = self.mechanism.gas
        jump = gas.properties(1000.+1.e-8, self.y)['h']-gas.properties(1000., self.y)['h']
        self.assertAlmostEqual(jump, -.139543, delta=1.e-5)
        # Starting above all 1000 K switches avoids this discontinuity.
        errors = []
        for tol in (1.e-5, 1.e-9):
            result = self.run_own(temperature=1100., rtol=tol,
                                 atol_species=tol*1.e-7, atol_temperature=tol*10.)
            errors.append(result.diagnostics['energy_relative_error'])
        self.assertLess(errors[1], errors[0]*.1)

    def test_negative_accepted_state_is_not_clipped(self):
        class InvalidBDF:
            def __init__(fake, fun, t0, state, bound, **kwargs):
                fake.t, fake.y, fake.status = t0, state.copy(), 'running'
            def step(fake):
                fake.t = 1.e-8
                fake.y[1] = -1.e-30
                fake.status = 'finished'
        with patch('scipy.integrate.BDF', InvalidBDF), self.assertRaisesRegex(MechanismError, 'not clipped'):
            self.run_own()

    def test_no_cantera_calls_during_integration(self):
        with patch.object(self.ct, 'Solution', side_effect=AssertionError('Cantera runtime call')):
            result = self.run_own(end_time=1.e-6)
        self.assertIsNone(result.ignition_delay)

    def test_invalid_inputs_and_step_limit(self):
        for kwargs in [dict(mode='unknown'), dict(rtol=0), dict(end_time=-1),
                       dict(max_step=0), dict(atol_species=0), dict(max_steps=True),
                       dict(mass_fractions=[.1]*9), dict(temperature=10),
                       dict(max_steps=1), dict(ignition_temperature_rise=float('nan'))]:
            with self.subTest(kwargs=kwargs), self.assertRaises(MechanismError):
                self.run_own(**kwargs)

    def test_stiff_analytic_and_inert(self):
        nasa = NASA('NASA7', (200, 4000), ((3.5, 0, 0, 0, 0, 0, 0),))
        gas = IdealGas(('A', 'B'), (.01, .01), (nasa, nasa))
        reaction = Reaction(Rate('Arrhenius', high=Arrhenius(math.log(1.e6), 0, 0)),
                            (1, 0), (0, 1), (1, 0), False)
        for reacting in (True, False):
            mechanism = SimpleNamespace(gas=gas, kinetics=Kinetics(gas, (reaction,) if reacting else ()),
                topology=SimpleNamespace(elements=('X',), compositions=((1,), (1,))))
            result = integrate(mechanism, temperature=1000., pressure=101325.,
                mass_fractions=[1., 0.], end_time=5.e-6, rtol=1.e-9, atol_species=1.e-16)
            self.assertAlmostEqual(result.mass_fractions[-1][0], math.exp(-5) if reacting else 1., delta=1.e-9)
            self.assertAlmostEqual(result.temperature[-1], 1000., places=8)
            self.assertIsNone(result.ignition_delay)

    def test_cli(self):
        root = Path(__file__).resolve().parents[1]
        command = [sys.executable, str(root/'tools/run_reference_reactor.py'), str(self.path),
                   '--temperature', '1000', '--mass-fractions', '{"H2":0.1,"O2":0.9}',
                   '--end-time', '1e-7']
        run = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        report = json.loads(run.stdout)
        self.assertEqual(report['result']['time'][-1], 1.e-7)
        self.assertEqual(len(report['source_sha256']), 64)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'reactor.json'
            saved = subprocess.run(command+['--output', str(output)], capture_output=True, text=True)
            self.assertEqual(saved.returncode, 0, saved.stderr)
            original = output.read_bytes()
            again = subprocess.run(command+['--output', str(output)], capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)
            self.assertEqual(output.read_bytes(), original)
            failed_output = Path(directory)/'failed.json'
            failed = subprocess.run(command+['--output', str(failed_output), '--max-steps', '1'],
                                    capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertFalse(failed_output.exists())
        bad = subprocess.run(command+['--mode', 'wrong'], capture_output=True, text=True)
        self.assertNotEqual(bad.returncode, 0)


if __name__ == '__main__':
    unittest.main()
