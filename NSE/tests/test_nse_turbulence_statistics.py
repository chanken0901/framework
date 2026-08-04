from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np


TOOLS_DIR = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from nse_turbulence_statistics import calculate_statistics  # noqa: E402


class TurbulenceStatisticsTests(unittest.TestCase):
    def test_single_solenoidal_wave(self) -> None:
        size = 16
        length = 2.0 * math.pi
        spacing = (length / size,) * 3
        z = np.arange(size, dtype=np.float64) * spacing[2]
        phase = 2.0 * z
        u = np.broadcast_to(np.sin(phase), (size, size, size)).copy()
        v = np.broadcast_to(np.cos(phase), (size, size, size)).copy()
        w = np.zeros_like(u)
        thermodynamics = {
            "mean_density": 1.0,
            "mean_pressure": 1.0 / 1.4,
            "mean_sound_speed": 1.0,
        }

        result = calculate_statistics(
            10,
            0.25,
            u,
            v,
            w,
            thermodynamics,
            spacing,
            reynolds_number=100.0,
        )

        self.assertAlmostEqual(result["rms_u_fluctuation"], 1.0 / math.sqrt(2.0), places=12)
        self.assertAlmostEqual(result["rms_v_fluctuation"], 1.0 / math.sqrt(2.0), places=12)
        self.assertAlmostEqual(result["rms_w_fluctuation"], 0.0, places=12)
        self.assertAlmostEqual(result["rms_velocity_fluctuation"], 1.0, places=12)
        self.assertAlmostEqual(result["turbulent_kinetic_energy"], 0.5, places=12)
        self.assertAlmostEqual(result["integral_length_scale"], 3.0 * math.pi / 8.0, places=12)
        self.assertAlmostEqual(result["dissipation_rate"], 0.04, places=12)
        self.assertAlmostEqual(result["taylor_microscale"], math.sqrt(1.25), places=12)
        self.assertAlmostEqual(
            result["kolmogorov_length_scale"], (0.01**3 / 0.04) ** 0.25, places=12
        )
        self.assertAlmostEqual(result["turbulent_mach_number"], 1.0, places=12)
        self.assertLess(result["spectral_energy_relative_error"], 1.0e-13)


if __name__ == "__main__":
    unittest.main()
