"""Validation failures must not silently pass the temperature regression."""
import unittest

import numpy as np

from check_watanabe_temperature import check_budget, check_samples, state_sample


class TemperatureValidationTests(unittest.TestCase):
    def sample(self):
        q = np.zeros((2, 2, 2, 5))
        q[..., 0] = 1
        q[..., 4] = 2.5
        return q

    def history(self):
        return [state_sample(self.sample(), t, 2) for t in np.linspace(0, .1, 11)]

    def log(self):
        return "\n".join(
            f"# forcing variance dissipation thermal_source .04 .012 -.012\n"
            f"# forcing {i} .3 0 .04 0 0 .012 0" for i in range(1, 31))

    def test_mass_and_volume_temperature_are_distinct(self):
        q = self.sample()
        q[0, ..., 0] = 2
        s = state_sample(q, 0, 2)
        self.assertAlmostEqual(s["mass_temperature"], 2/3)
        self.assertAlmostEqual(s["volume_temperature"], .75)

    def test_invalid_states(self):
        for component, value in [(0, 0), (0, -1), (4, -1), (4, np.nan), (1, 1e308)]:
            with self.subTest(component=component, value=value):
                q = self.sample()
                q[0, 0, 0, component] = value
                with self.assertRaises((ValueError, FloatingPointError)):
                    state_sample(q, 0, 2)
        with self.assertRaises(ValueError):
            state_sample(self.sample(), np.nan, 2)

    def test_ghost_crop(self):
        q = np.full((8, 8, 8, 5), np.nan)
        q[3:-3, 3:-3, 3:-3] = self.sample()
        self.assertEqual(state_sample(q, 0, 2)["mass_temperature"], 1)

    def test_output_coverage(self):
        history = self.history()
        self.assertEqual(check_samples(history, .1)["mass_temperature"], 0)
        for bad in [history[:-1], history[:4]+history[5:], history+[history[-1]],
                    history[:4]+[history[3]]+history[5:]]:
            with self.assertRaises(ValueError):
                check_samples(bad, .1)

    def test_budget_and_missing_heat_work(self):
        self.assertLess(check_budget(self.log(), 10, .01, self.history())[
            "max_energy_budget_residual"], 1e-15)
        for rate in [.012, -.012]:
            history = self.history()
            for s in history:
                s["energy"] += rate*s["time"]
            with self.assertRaisesRegex(ValueError, "budget mismatch"):
                check_budget(self.log(), 10, .01, history)

    def test_bad_diagnostics(self):
        for log in ["", self.log().replace("-.012", "0"),
                    self.log().replace(".3 0", ".2 0"),
                    self.log().replace(".04 .012", "nan .012"),
                    self.log().replace("# forcing 2 ", "# forcing 1 ")]:
            with self.assertRaises(ValueError):
                check_budget(log, 10, .01, self.history())

    def test_ssprk_stage_weights(self):
        history = self.history()
        lines = []
        rates = []
        for i in range(30):
            variance = [.04, .05, .06][i % 3]
            coefficient = (variance/.04)**(-4)*.012/.04
            rates.append(coefficient*.04-.012)
            lines.extend([
                f"# forcing {i+1} {coefficient:.16e} 0 .04 0 0 .012 0",
                f"# forcing variance dissipation thermal_source {variance} .012 -.012",
            ])
        increment = .01*(rates[0]+rates[1]+4*rates[2])/6
        for i, sample in enumerate(history):
            sample["energy"] += i*increment
        budget = check_budget("\n".join(lines), 10, .01, history)
        self.assertLess(budget["max_energy_budget_residual"], 1e-14)


if __name__ == "__main__":
    unittest.main()
