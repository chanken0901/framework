from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
FRAMEWORK_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from case_input import (  # noqa: E402
    CaseInputError,
    NSE_BOUNDARY_FACES,
    _validate_solver_selection,
    derive_nse_hit_transport,
    render_nse,
    render_nse_multicomponent,
)
from yaml_support import load_yaml  # noqa: E402


MANIFEST = (
    FRAMEWORK_ROOT / "SolverLibrary" / "GPE" / "gp3d" / "solver_manifest.yaml"
)
NSE_MANIFEST = FRAMEWORK_ROOT / "SolverLibrary" / "NSE" / "solver_manifest.yaml"


class CaseInputProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(MANIFEST)

    def test_rejects_profile_that_differs_from_case(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_cuda": False,
                "mpi_processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "does not match"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_rejects_multiple_mpi_processes_for_single_gpu(self) -> None:
        case = {
            "solver": {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_cuda": True,
                "mpi_processes": 8,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "must be 1"):
            _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_accepts_matching_cpu_mpi_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_cuda": False,
                "mpi_processes": 8,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_profile_derived_fields_may_be_omitted(self) -> None:
        case = {
            "solver": {
                "use_openmp": True,
                "mpi_processes": 8,
                "omp_threads": 2,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_accepts_runtime_openmp_for_hybrid_capable_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cpu_mpi_fftw",
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
                "mpi_processes": 4,
                "omp_threads": 4,
            }
        }

        _validate_solver_selection(case, self.manifest, "cpu_mpi_fftw")

    def test_rejects_runtime_openmp_for_non_openmp_profile(self) -> None:
        case = {
            "solver": {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": True,
                "use_cuda": True,
                "mpi_processes": 1,
                "omp_threads": 2,
            }
        }

        with self.assertRaisesRegex(CaseInputError, "without OpenMP"):
            _validate_solver_selection(case, self.manifest, "cuda_single")


class MulticomponentFoundationInputTests(unittest.TestCase):
    def test_stage_eight_pencil_namelist(self) -> None:
        case = self.reactive_boundary_case()
        case["solver"] = {"mpi_processes": 4, "decomposition": "pencil", "process_grid": [2, 2]}
        text = render_nse_multicomponent(case, "cpu_mpi_reactive_pencil")
        self.assertIn("&multicomponent_parallel", text)
        self.assertIn('decomposition = "pencil"', text)
        self.assertIn("process_grid = 2, 2", text)
        self.assertIn("boundary_face_types", text)

    def test_stage_eight_rejects_invalid_pencil_grid(self) -> None:
        for grid in ([3, 2], [-1, 0], [True, 4], [2], [3, 0]):
            with self.subTest(grid=grid):
                case = self.reactive_boundary_case()
                case["solver"] = {"mpi_processes": 4, "process_grid": grid}
                with self.assertRaises(CaseInputError):
                    render_nse_multicomponent(case, "cpu_mpi_reactive_pencil")

    def test_stage_eight_rejects_slab(self) -> None:
        case = self.reactive_boundary_case()
        case["solver"] = {"decomposition": "slab"}
        with self.assertRaisesRegex(CaseInputError, "decomposition=pencil"):
            render_nse_multicomponent(case, "cpu_mpi_reactive_pencil")

    def test_stage_eight_openmp_retains_physics_input(self) -> None:
        case = self.reactive_boundary_case()
        self.assertEqual(
            render_nse_multicomponent(case, "cpu_openmp_reactive"),
            render_nse_multicomponent(case, "cpu_serial_reactive_boundaries"),
        )

    @staticmethod
    def case() -> dict:
        return {
            "physics": {
                "model": "nse_multicomponent",
                "multicomponent": {"species": ["mixture"]},
            },
            "thermodynamics": {"model": "calorically_perfect"},
            "transport": {"model": "none"},
            "chemistry": {"model": "none"},
        }

    @classmethod
    def euler_case(cls) -> dict:
        case = cls.case()
        case["physics"]["multicomponent"] = {
            "mode": "inviscid_euler",
            "species": ["species_a", "species_b"],
        }
        case["thermodynamics"]["gamma"] = 1.4
        case["grid"] = {
            "nx": 32,
            "ny": 4,
            "nz": 4,
            "x_min": 0.0,
            "x_max": 1.0,
            "y_min": 0.0,
            "y_max": 1.0,
            "z_min": 0.0,
            "z_max": 1.0,
        }
        case["flow"] = {
            "type": "multispecies_sod",
            "multispecies_sod": {
                "initial_condition": "multispecies_sod_x",
                "interface_location": 0.5,
                "left": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 1.0,
                    "mass_fractions": [0.8, 0.2],
                },
                "right": {
                    "density": 0.125,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 0.1,
                    "mass_fractions": [0.2, 0.8],
                },
            },
        }
        case["time"] = {"cfl": 0.35, "dt": 0.0, "nsteps": 10}
        case["numerics"] = {
            "convective_scheme": "rusanov1",
            "boundary_condition": "periodic",
            "time_integration": "ssprk3",
        }
        case["output"] = {
            "write_final": True,
            "filename": "multicomponent_euler_final.csv",
        }
        return case

    @classmethod
    def thermally_perfect_case(cls) -> dict:
        case = cls.euler_case()
        case["physics"]["multicomponent"] = {
            "mode": "thermally_perfect_euler",
            "species": ["N2", "O2"],
        }
        case["thermodynamics"] = {
            "model": "thermally_perfect",
            "universal_gas_constant": 8314.46261815324,
            "temperature_min": 200.0,
            "temperature_max": 6000.0,
            "species_data": {
                "N2": {
                    "molecular_weight": 28.0134,
                    "temperature_midpoint": 1000.0,
                    "nasa7_low": [
                        3.53100528,
                        -1.23660987e-4,
                        -5.02999433e-7,
                        2.43530612e-9,
                        -1.40881235e-12,
                        -1046.97628,
                        2.96747468,
                    ],
                    "nasa7_high": [
                        2.95257626,
                        1.39690040e-3,
                        -4.92631603e-7,
                        7.86010367e-11,
                        -4.60755321e-15,
                        -923.948645,
                        5.87188762,
                    ],
                },
                "O2": {
                    "molecular_weight": 31.9988,
                    "temperature_midpoint": 1000.0,
                    "nasa7_low": [
                        3.78245636,
                        -2.99673416e-3,
                        9.84730201e-6,
                        -9.68129509e-9,
                        3.24372837e-12,
                        -1063.94356,
                        3.65767573,
                    ],
                    "nasa7_high": [
                        3.28253784,
                        1.48308754e-3,
                        -7.57966669e-7,
                        2.09470555e-10,
                        -2.16717794e-14,
                        -1088.45772,
                        5.45323129,
                    ],
                },
            },
        }
        case["flow"]["multispecies_sod"]["left"] = {
            "density": 1.17197031944841,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 101325.0,
            "mass_fractions": [0.767, 0.233],
        }
        case["flow"]["multispecies_sod"]["right"] = {
            "density": 0.126389266436632,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 10132.5,
            "mass_fractions": [0.2, 0.8],
        }
        case["output"]["filename"] = (
            "multicomponent_thermally_perfect_final.csv"
        )
        return case

    @classmethod
    def viscous_case(cls) -> dict:
        case = cls.thermally_perfect_case()
        case["physics"]["multicomponent"]["mode"] = (
            "viscous_navier_stokes"
        )
        case["transport"] = {
            "model": "mixture_averaged",
            "reference_dynamic_viscosity": 1.8e-5,
            "prandtl_number": 0.72,
            "species_data": {
                "N2": {"diffusivity": 2.0e-5},
                "O2": {"diffusivity": 2.1e-5},
            },
        }
        case["flow"] = {
            "type": "periodic_species_wave",
            "periodic_species_wave": {
                "initial_condition": "periodic_species_wave_x",
                "density": 1.0,
                "temperature": 300.0,
                "velocity": [0.0, 0.0, 0.0],
                "mean_mass_fractions": [0.5, 0.5],
                "positive_species": "N2",
                "negative_species": "O2",
                "amplitude": 0.1,
                "wavenumber": 1,
            },
        }
        case["time"] = {
            "cfl": 0.2,
            "diffusion_cfl": 0.4,
            "dt": 0.0,
            "nsteps": 10,
        }
        case["output"]["filename"] = "multicomponent_viscous_final.csv"
        return case

    @classmethod
    def reactor_case(cls) -> dict:
        case = cls.case()
        case["physics"]["multicomponent"] = {
            "mode": "homogeneous_reactor",
            "species": ["fuel", "oxidizer", "product"],
        }
        species_data = {}
        for name, formation_coefficient in (
            ("fuel", 0.0),
            ("oxidizer", 0.0),
            ("product", -5000.0),
        ):
            coefficients = [3.5, 0.0, 0.0, 0.0, 0.0, formation_coefficient, 0.0]
            species_data[name] = {
                "molecular_weight": 28.0,
                "temperature_midpoint": 1000.0,
                "nasa7_low": coefficients,
                "nasa7_high": coefficients,
            }
        case["thermodynamics"] = {
            "model": "thermally_perfect",
            "universal_gas_constant": 8314.46261815324,
            "temperature_min": 200.0,
            "temperature_max": 6000.0,
            "species_data": species_data,
        }
        case["chemistry"] = {
            "model": "one_step_arrhenius",
            "reaction": {
                "reactants": {"fuel": 1.0, "oxidizer": 1.0},
                "products": {"product": 2.0},
                "pre_exponential_factor": 1000.0,
                "temperature_exponent": 0.0,
                "activation_temperature": 2000.0,
            },
        }
        case["flow"] = {
            "type": "homogeneous_reactor",
            "homogeneous_reactor": {
                "density": 1.0,
                "temperature": 1200.0,
                "mass_fractions": {
                    "fuel": 0.5,
                    "oxidizer": 0.5,
                    "product": 0.0,
                },
            },
        }
        case["time"] = {
            "dt": 0.0,
            "maximum_dt": 1.0e-3,
            "chemistry_cfl": 0.1,
            "nsteps": 100,
        }
        case["output"] = {
            "write_history": True,
            "output_every": 5,
            "filename": "homogeneous_reactor.csv",
        }
        return case

    @classmethod
    def reactive_case(cls) -> dict:
        case = cls.reactor_case()
        case["physics"]["multicomponent"]["mode"] = (
            "reactive_navier_stokes"
        )
        case["transport"] = {
            "model": "mixture_averaged",
            "reference_dynamic_viscosity": 1.8e-5,
            "prandtl_number": 0.72,
            "species_data": {
                "fuel": {"diffusivity": 2.0e-5},
                "oxidizer": {"diffusivity": 2.0e-5},
                "product": {"diffusivity": 2.0e-5},
            },
        }
        case["grid"] = {
            "nx": 24,
            "ny": 4,
            "nz": 4,
            "x_min": 0.0,
            "x_max": 1.0,
            "y_min": 0.0,
            "y_max": 1.0,
            "z_min": 0.0,
            "z_max": 1.0,
        }
        case["flow"] = {
            "type": "periodic_species_wave",
            "periodic_species_wave": {
                "initial_condition": "periodic_species_wave_x",
                "density": 1.0,
                "temperature": 1200.0,
                "velocity": [0.0, 0.0, 0.0],
                "mean_mass_fractions": [0.45, 0.45, 0.10],
                "positive_species": "fuel",
                "negative_species": "product",
                "amplitude": 0.05,
                "wavenumber": 1,
            },
        }
        case["time"] = {
            "cfl": 0.2,
            "diffusion_cfl": 0.4,
            "chemistry_cfl": 0.1,
            "maximum_chemistry_substeps": 10000,
            "dt": 0.0,
            "nsteps": 20,
        }
        case["numerics"] = {
            "convective_scheme": "rusanov1",
            "boundary_condition": "periodic",
            "time_integration": "ssprk3",
            "coupling_scheme": "strang",
            "chemistry_time_integration": "ssprk3_subcycled",
        }
        case["output"] = {
            "write_final": True,
            "filename": "multicomponent_reactive_final.csv",
        }
        return case

    @classmethod
    def reactive_boundary_case(cls) -> dict:
        case = cls.reactive_case()
        case["flow"] = {
            "type": "reactive_shock_tube",
            "reactive_shock_tube": {
                "initial_condition": "reactive_shock_tube_x",
                "interface_location": 0.35,
                "left": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 356334.11220656743,
                    "mass_fractions": {
                        "fuel": 0.45,
                        "oxidizer": 0.45,
                        "product": 0.10,
                    },
                },
                "right": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 267250.5841549256,
                    "mass_fractions": {
                        "fuel": 0.49,
                        "oxidizer": 0.49,
                        "product": 0.02,
                    },
                },
            },
        }
        case["numerics"].pop("boundary_condition")
        case["boundary"] = {
            "faces": {
                "x_min": {
                    "type": "dirichlet",
                    "reference_state": "driver",
                },
                "x_max": {
                    "type": "non_reflecting",
                    "reference_state": "far_field",
                },
                "y_min": {"type": "periodic"},
                "y_max": {"type": "periodic"},
                "z_min": {"type": "reflective"},
                "z_max": {"type": "reflective"},
            },
            "reference_states": {
                "driver": case["flow"]["reactive_shock_tube"]["left"],
                "far_field": case["flow"]["reactive_shock_tube"]["right"],
            },
            "non_reflecting": {
                "formulation": "characteristic_relaxation",
                "relaxation_strength": 0.1,
                "length_scale": "auto",
            },
        }
        case["output"].update(
            {
                "write_snapshots": True,
                "write_history": True,
                "output_every": 5,
                "snapshot_prefix": "reactive_shock_tube",
                "history_filename": "reactive_shock_tube_history.csv",
            }
        )
        return case

    def test_renders_one_species_stage_zero_contract(self) -> None:
        text = render_nse_multicomponent(self.case())

        self.assertIn("&multicomponent", text)
        self.assertIn("nspecies = 1", text)
        self.assertIn('species_names = "mixture"', text)
        self.assertIn('simulation_mode = "foundation"', text)
        self.assertIn('thermodynamics_model = "calorically_perfect"', text)
        self.assertIn('transport_model = "none"', text)
        self.assertIn('chemistry_model = "none"', text)

    def test_renders_multiple_species_without_changing_the_contract(self) -> None:
        case = self.case()
        case["physics"]["multicomponent"]["species"] = ["fuel", "oxidizer"]

        text = render_nse_multicomponent(case)

        self.assertIn("nspecies = 2", text)
        self.assertIn('species_names = "fuel", "oxidizer"', text)

    def test_rejects_unimplemented_reactive_provider(self) -> None:
        case = self.case()
        case["chemistry"]["model"] = "finite_rate"

        with self.assertRaisesRegex(CaseInputError, "multicomponent.*chemistry"):
            render_nse_multicomponent(case)

    def test_rejects_duplicate_species(self) -> None:
        case = self.case()
        case["physics"]["multicomponent"]["species"] = ["N2", "N2"]

        with self.assertRaisesRegex(CaseInputError, "unique"):
            render_nse_multicomponent(case)

    def test_renders_stage_one_passive_scalar_contract(self) -> None:
        case = self.case()
        case["physics"]["multicomponent"] = {
            "mode": "passive_scalar",
            "species": ["tracer", "carrier"],
        }
        case.update(
            {
                "grid": {
                    "nx": 16,
                    "ny": 8,
                    "nz": 4,
                    "x_min": 0.0,
                    "x_max": 1.0,
                    "y_min": 0.0,
                    "y_max": 1.0,
                    "z_min": 0.0,
                    "z_max": 1.0,
                },
                "flow": {
                    "type": "passive_scalar_advection",
                    "velocity": [1.0, -0.25, 0.0],
                    "passive_scalar": {
                        "initial_condition": "gaussian",
                        "tracer_center": [0.25, 0.5, 0.5],
                    },
                },
                "time": {"cfl": 0.45, "dt": 0.0, "nsteps": 20},
                "numerics": {
                    "convective_scheme": "upwind1",
                    "boundary_condition": "periodic",
                    "time_integration": "ssprk3",
                },
                "output": {
                    "write_final": True,
                    "filename": "passive_scalar_final.csv",
                },
            }
        )

        text = render_nse_multicomponent(case)

        self.assertIn('simulation_mode = "passive_scalar"', text)
        self.assertIn("&passive_scalar", text)
        self.assertIn("velocity = 1, -0.25, 0", text)
        self.assertIn('advection_scheme = "upwind1"', text)
        self.assertIn('boundary_condition = "periodic"', text)
        self.assertIn('time_integrator = "ssprk3"', text)

    def test_passive_scalar_requires_two_species(self) -> None:
        case = self.case()
        case["physics"]["multicomponent"]["mode"] = "passive_scalar"

        with self.assertRaisesRegex(CaseInputError, "tracer and carrier"):
            render_nse_multicomponent(case)

    def test_profile_must_match_multicomponent_mode(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "requires.*passive_scalar"):
            render_nse_multicomponent(
                self.case(), "cpu_serial_passive_scalar"
            )

    def test_renders_stage_two_multicomponent_euler_contract(self) -> None:
        text = render_nse_multicomponent(
            self.euler_case(), "cpu_serial_inviscid"
        )

        self.assertIn('simulation_mode = "inviscid_euler"', text)
        self.assertIn("&multicomponent_euler", text)
        self.assertIn("gamma = 1.3999999999999999", text)
        self.assertIn("left_mass_fractions = 0.80000000000000004", text)
        self.assertIn('riemann_solver = "rusanov1"', text)
        self.assertIn('boundary_condition = "periodic"', text)
        self.assertIn('time_integrator = "ssprk3"', text)

    def test_euler_rejects_mass_fractions_that_do_not_sum_to_one(self) -> None:
        case = self.euler_case()
        case["flow"]["multispecies_sod"]["left"]["mass_fractions"] = [
            0.8,
            0.3,
        ]

        with self.assertRaisesRegex(CaseInputError, "must sum to one"):
            render_nse_multicomponent(case)

    def test_euler_profile_rejects_foundation_mode(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "requires.*inviscid_euler"):
            render_nse_multicomponent(self.case(), "cpu_serial_inviscid")

    def test_renders_stage_three_thermally_perfect_contract(self) -> None:
        text = render_nse_multicomponent(
            self.thermally_perfect_case(),
            "cpu_serial_thermally_perfect",
        )

        self.assertIn('simulation_mode = "thermally_perfect_euler"', text)
        self.assertIn('thermodynamics_model = "thermally_perfect"', text)
        self.assertIn("&thermally_perfect", text)
        self.assertIn('thermo_species_names = "N2", "O2"', text)
        self.assertIn("molecular_weights = 28.0134", text)
        self.assertIn("nasa_low_coefficients = 3.53100528", text)
        self.assertIn("&multicomponent_euler", text)

    def test_stage_three_requires_exact_species_property_keys(self) -> None:
        case = self.thermally_perfect_case()
        del case["thermodynamics"]["species_data"]["O2"]

        with self.assertRaisesRegex(CaseInputError, "must exactly match"):
            render_nse_multicomponent(case)

    def test_stage_three_rejects_initial_temperature_outside_range(self) -> None:
        case = self.thermally_perfect_case()
        case["flow"]["multispecies_sod"]["left"]["density"] = 0.01

        with self.assertRaisesRegex(CaseInputError, "outside.*temperature range"):
            render_nse_multicomponent(case)

    def test_stage_three_profile_rejects_stage_two_mode(self) -> None:
        with self.assertRaisesRegex(
            CaseInputError, "requires.*thermally_perfect_euler"
        ):
            render_nse_multicomponent(
                self.euler_case(), "cpu_serial_thermally_perfect"
            )

    def test_renders_stage_four_transport_contract(self) -> None:
        text = render_nse_multicomponent(
            self.viscous_case(), "cpu_serial_viscous"
        )

        self.assertIn('simulation_mode = "viscous_navier_stokes"', text)
        self.assertIn('transport_model = "mixture_averaged"', text)
        self.assertIn("&mixture_averaged_transport", text)
        self.assertIn('transport_species_names = "N2", "O2"', text)
        self.assertIn("species_diffusivities = 2.0000000000000002e-05", text)
        self.assertIn('initial_condition = "periodic_species_wave_x"', text)
        self.assertIn("wave_positive_species = 1", text)
        self.assertIn("wave_negative_species = 2", text)
        self.assertIn("diffusion_cfl = 0.40000000000000002", text)

    def test_stage_four_transport_species_keys_must_match(self) -> None:
        case = self.viscous_case()
        del case["transport"]["species_data"]["O2"]

        with self.assertRaisesRegex(CaseInputError, "must exactly match"):
            render_nse_multicomponent(case, "cpu_serial_viscous")

    def test_stage_four_rejects_nonpositive_diffusivity(self) -> None:
        case = self.viscous_case()
        case["transport"]["species_data"]["O2"]["diffusivity"] = 0.0

        with self.assertRaisesRegex(CaseInputError, "must be positive"):
            render_nse_multicomponent(case, "cpu_serial_viscous")

    def test_stage_four_rejects_wave_that_breaks_positivity(self) -> None:
        case = self.viscous_case()
        case["flow"]["periodic_species_wave"]["amplitude"] = 0.5

        with self.assertRaisesRegex(CaseInputError, "violates.*positivity"):
            render_nse_multicomponent(case, "cpu_serial_viscous")

    def test_renders_stage_five_homogeneous_reactor_contract(self) -> None:
        text = render_nse_multicomponent(
            self.reactor_case(), "cpu_serial_reactor"
        )

        self.assertIn('simulation_mode = "homogeneous_reactor"', text)
        self.assertIn('chemistry_model = "one_step_arrhenius"', text)
        self.assertIn("&one_step_arrhenius", text)
        self.assertIn("reactant_stoich = 1, 1, 0", text)
        self.assertIn("product_stoich = 0, 0, 2", text)
        self.assertIn("reaction_orders = 1, 1, 0", text)
        self.assertIn("&homogeneous_reactor", text)
        self.assertIn("initial_mass_fractions = 0.5, 0.5, 0", text)
        self.assertIn("chemistry_cfl = 0.1", text)
        self.assertIn('output_file = "homogeneous_reactor.csv"', text)

    def test_stage_five_rejects_nonconservative_stoichiometry(self) -> None:
        case = self.reactor_case()
        case["chemistry"]["reaction"]["products"]["product"] = 1.0

        with self.assertRaisesRegex(CaseInputError, "does not conserve mass"):
            render_nse_multicomponent(case, "cpu_serial_reactor")

    def test_stage_five_mass_fractions_are_named_and_complete(self) -> None:
        case = self.reactor_case()
        del case["flow"]["homogeneous_reactor"]["mass_fractions"]["product"]

        with self.assertRaisesRegex(CaseInputError, "keys must exactly match"):
            render_nse_multicomponent(case, "cpu_serial_reactor")

    def test_renders_stage_six_reactive_flow_contract(self) -> None:
        text = render_nse_multicomponent(
            self.reactive_case(), "cpu_serial_reactive"
        )

        self.assertIn('simulation_mode = "reactive_navier_stokes"', text)
        self.assertIn('transport_model = "mixture_averaged"', text)
        self.assertIn('chemistry_model = "one_step_arrhenius"', text)
        self.assertIn("&thermally_perfect", text)
        self.assertIn("&mixture_averaged_transport", text)
        self.assertIn("&one_step_arrhenius", text)
        self.assertIn("&multicomponent_euler", text)
        self.assertIn("&reactive_navier_stokes", text)
        self.assertIn('splitting_scheme = "strang"', text)
        self.assertIn(
            'chemistry_integrator = "ssprk3_subcycled"', text
        )
        self.assertIn("maximum_chemistry_substeps = 10000", text)
        self.assertIn(
            'output_file = "multicomponent_reactive_final.csv"', text
        )

    def test_stage_six_rejects_unsupported_coupling_scheme(self) -> None:
        case = self.reactive_case()
        case["numerics"]["coupling_scheme"] = "lie"

        with self.assertRaisesRegex(CaseInputError, "coupling_scheme"):
            render_nse_multicomponent(case, "cpu_serial_reactive")

    def test_stage_six_requires_reactive_profile_mode(self) -> None:
        with self.assertRaisesRegex(
            CaseInputError, "requires.*reactive_navier_stokes"
        ):
            render_nse_multicomponent(
                self.viscous_case(), "cpu_serial_reactive"
            )

    def test_renders_stage_seven_boundary_and_output_contract(self) -> None:
        text = render_nse_multicomponent(
            self.reactive_boundary_case(),
            "cpu_serial_reactive_boundaries",
        )

        self.assertIn('initial_condition = "reactive_shock_tube_x"', text)
        self.assertIn('boundary_condition = "face_specific"', text)
        self.assertIn(
            'boundary_face_types = "dirichlet", "non_reflecting", '
            '"periodic", "periodic", "reflective", "reflective"',
            text,
        )
        self.assertIn(
            "boundary_reference_mass_fractions = 0.45000000000000001",
            text,
        )
        self.assertIn("write_snapshots = .true.", text)
        self.assertIn("write_history = .true.", text)
        self.assertIn("output_every = 5", text)
        self.assertIn('history_file = "reactive_shock_tube_history.csv"', text)

    def test_stage_seven_rejects_unpaired_periodic_faces(self) -> None:
        case = self.reactive_boundary_case()
        case["boundary"]["faces"]["y_max"] = {"type": "reflective"}

        with self.assertRaisesRegex(CaseInputError, "periodic y.*paired"):
            render_nse_multicomponent(
                case, "cpu_serial_reactive_boundaries"
            )

    def test_stage_seven_requires_boundary_reference_composition(self) -> None:
        case = self.reactive_boundary_case()
        case["boundary"]["reference_states"]["far_field"] = dict(
            case["boundary"]["reference_states"]["far_field"]
        )
        del case["boundary"]["reference_states"]["far_field"][
            "mass_fractions"
        ]

        with self.assertRaisesRegex(CaseInputError, "missing: mass_fractions"):
            render_nse_multicomponent(
                case, "cpu_serial_reactive_boundaries"
            )

    def test_stage_six_profile_rejects_stage_seven_features(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "Stage-6 periodic profile"):
            render_nse_multicomponent(
                self.reactive_boundary_case(), "cpu_serial_reactive"
            )

    def test_passive_scalar_rejects_unimplemented_scheme(self) -> None:
        case = self.case()
        case["physics"]["multicomponent"] = {
            "mode": "passive_scalar",
            "species": ["tracer", "carrier"],
        }
        case["grid"] = {
            "nx": 4,
            "ny": 4,
            "nz": 4,
            "x_min": 0.0,
            "x_max": 1.0,
            "y_min": 0.0,
            "y_max": 1.0,
            "z_min": 0.0,
            "z_max": 1.0,
        }
        case["flow"] = {
            "type": "passive_scalar_advection",
            "velocity": [1.0, 0.0, 0.0],
            "passive_scalar": {},
        }
        case["time"] = {"nsteps": 1}
        case["numerics"] = {"convective_scheme": "weno5z_roe"}

        with self.assertRaisesRegex(CaseInputError, "requires.*upwind1"):
            render_nse_multicomponent(case)


class NseCaseInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = load_yaml(NSE_MANIFEST)

    @staticmethod
    def case(mpi_processes: int = 4) -> dict:
        return {
            "case_id": "case0001",
            "physics": {
                "model": "nse",
                "nse": {
                    "nv": 5,
                    "gamma": 1.4,
                    "rho0": 1.0,
                    "mach_number": 0.5,
                    "reynolds_number": 0.0,
                    "prandtl_number": 0.72,
                },
            },
            "flow": {"type": "tgv"},
            "grid": {
                "nx": 8,
                "ny": 8,
                "nz": 8,
                "nghost": 3,
                "x_min": 0.0,
                "x_max": 6.283185307179586,
                "y_min": 0.0,
                "y_max": 6.283185307179586,
                "z_min": 0.0,
                "z_max": 6.283185307179586,
            },
            "time": {
                "dt": 1.0e-4,
                "t_max": 1.0e-4,
                "nsteps": 1,
                "cfl": 0.5,
                "use_fixed_dt": False,
                "output_frequency": 1,
            },
            "solver": {
                "profile": "cpu_mpi",
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
                "mpi_processes": mpi_processes,
                "omp_threads": 2,
            },
            "numerics": {
                "convective_scheme": "KEEP6",
                "viscous_scheme": "none",
                "boundary_condition": "periodic",
                "time_integration": "RK3",
            },
            "output": {},
        }

    @staticmethod
    def x_non_reflecting_boundary() -> dict:
        return {
            "faces": {
                "x_min": {
                    "type": "non_reflecting",
                    "reference_state": "far_field",
                },
                "x_max": {
                    "type": "non_reflecting",
                    "reference_state": "far_field",
                },
                "y_min": {"type": "periodic"},
                "y_max": {"type": "periodic"},
                "z_min": {"type": "periodic"},
                "z_max": {"type": "periodic"},
            },
            "reference_states": {
                "far_field": {
                    "density": 1.0,
                    "velocity": [0.5, 0.0, 0.0],
                    "pressure": 1.0 / 1.4,
                }
            },
            "non_reflecting": {
                "formulation": "characteristic_relaxation",
                "relaxation_strength": 0.1,
                "length_scale": "auto",
            },
        }

    def shock_turbulence_case(self) -> dict:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        dx = (case["grid"]["x_max"] - case["grid"]["x_min"]) / case[
            "grid"
        ]["nx"]
        case["flow"] = {
            "type": "shock_turbulence_interaction",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "mode": "embed",
                "x_start": case["grid"]["x_min"] + 4 * dx,
                "blend_cells": 0,
            },
            "planar_shock": {
                "position": case["grid"]["x_min"] + 2 * dx,
                "propagation_direction": "positive_x",
                "mach_number": 1.5,
                "upstream": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 1.0 / 1.4,
                },
            },
        }
        case["boundary"] = {
            "faces": {
                "x_min": {
                    "type": "dirichlet",
                    "reference_state": "shock_post",
                },
                "x_max": {
                    "type": "non_reflecting",
                    "reference_state": "shock_pre",
                },
                "y_min": {"type": "periodic"},
                "y_max": {"type": "periodic"},
                "z_min": {"type": "periodic"},
                "z_max": {"type": "periodic"},
            },
            "reference_states": {
                "shock_pre": {"source": "planar_shock.upstream"},
                "shock_post": {"source": "planar_shock.downstream"},
            },
        }
        return case

    def shock_tube_turbulence_case(self) -> dict:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        dx = (case["grid"]["x_max"] - case["grid"]["x_min"]) / case[
            "grid"
        ]["nx"]
        case["flow"] = {
            "type": "shock_tube_turbulence_interaction",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "mode": "embed",
                "x_start": case["grid"]["x_min"] + 4 * dx,
                "blend_cells": 0,
            },
            "shock_tube": {
                "diaphragm_position": case["grid"]["x_min"] + 2 * dx,
                "driver": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 5.0 / 1.4,
                },
                "driven": {
                    "density": 1.0,
                    "velocity": [0.0, 0.0, 0.0],
                    "pressure": 1.0 / 1.4,
                },
            },
        }
        case["boundary"] = {
            "faces": {
                "x_min": {"type": "reflective"},
                "x_max": {
                    "type": "non_reflecting",
                    "reference_state": "driven",
                },
                "y_min": {"type": "periodic"},
                "y_max": {"type": "periodic"},
                "z_min": {"type": "periodic"},
                "z_max": {"type": "periodic"},
            },
            "reference_states": {
                "driven": {"source": "shock_tube.driven"},
            },
        }
        return case

    def test_renders_current_modular_nse_input(self) -> None:
        text = render_nse(self.case(), self.manifest, "cpu_mpi")

        self.assertIn('initial_condition = "taylor_green"', text)
        self.assertIn('backend = "cpu_mpi"', text)
        self.assertIn('convective_scheme = "keep6"', text)
        self.assertNotIn("hybrid_smooth_scheme", text)
        self.assertNotIn("hybrid_shock_scheme", text)
        self.assertNotIn("convective_order", text)
        self.assertIn('viscous_scheme = "none"', text)
        self.assertIn('boundary_condition = "periodic"', text)
        self.assertIn('time_integrator = "ssprk3"', text)

    def test_renders_x_non_reflecting_yz_periodic_boundaries(self) -> None:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        case["boundary"] = self.x_non_reflecting_boundary()

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertNotIn("boundary_condition =", text)
        self.assertIn('boundary_x_min = "non_reflecting"', text)
        self.assertIn('boundary_x_max = "non_reflecting"', text)
        self.assertIn('boundary_y_min = "periodic"', text)
        self.assertIn('boundary_z_max = "periodic"', text)
        self.assertIn("boundary_x_min_reference_rho = 1", text)
        self.assertIn("boundary_x_min_reference_u = 0.5", text)
        self.assertIn("boundary_relaxation_strength = 0.10000000000000001", text)
        self.assertIn("boundary_length_scale = -1", text)

    def test_rejects_unpaired_periodic_boundary(self) -> None:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        case["boundary"] = self.x_non_reflecting_boundary()
        case["boundary"]["faces"]["y_max"] = {
            "type": "non_reflecting",
            "reference_state": "far_field",
        }

        with self.assertRaisesRegex(CaseInputError, "periodic y boundaries"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_invalid_unreferenced_boundary_state(self) -> None:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        boundary = self.x_non_reflecting_boundary()
        boundary["reference_states"]["unused"] = {
            "density": -1.0,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 1.0,
        }
        case["boundary"] = boundary

        with self.assertRaisesRegex(CaseInputError, "unused.density"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_new_and_legacy_boundary_settings_together(self) -> None:
        case = self.case()
        case["boundary"] = self.x_non_reflecting_boundary()

        with self.assertRaisesRegex(CaseInputError, "cannot be specified together"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_non_reflecting_boundary_for_cuda(self) -> None:
        case = self.case(mpi_processes=1)
        case["numerics"].pop("boundary_condition")
        case["boundary"] = self.x_non_reflecting_boundary()
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )

        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('boundary_x_min = "non_reflecting"', text)
        self.assertIn('boundary_x_max = "non_reflecting"', text)
        self.assertIn('boundary_y_min = "periodic"', text)

    def test_accepts_new_all_periodic_boundary_for_cuda(self) -> None:
        case = self.case(mpi_processes=1)
        case["numerics"].pop("boundary_condition")
        boundary = self.x_non_reflecting_boundary()
        boundary["faces"] = {
            face: {"type": "periodic"} for face in NSE_BOUNDARY_FACES
        }
        boundary["reference_states"] = {}
        case["boundary"] = boundary
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )

        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('boundary_x_min = "periodic"', text)
        self.assertIn('boundary_z_max = "periodic"', text)
        self.assertNotIn("boundary_condition =", text)

    def test_renders_reflective_boundary_for_cpu_and_cuda(self) -> None:
        for profile, use_mpi, use_openmp, use_cuda in (
            ("cpu_mpi", True, True, False),
            ("cuda_single", False, False, True),
        ):
            with self.subTest(profile=profile):
                case = self.case(mpi_processes=4 if use_mpi else 1)
                case["numerics"].pop("boundary_condition")
                case["boundary"] = {
                    "faces": {
                        face: {"type": "reflective"}
                        for face in NSE_BOUNDARY_FACES
                    },
                    "reference_states": {},
                }
                case["solver"].update(
                    {
                        "profile": profile,
                        "use_mpi": use_mpi,
                        "use_openmp": use_openmp,
                        "use_cuda": use_cuda,
                        "omp_threads": 2 if use_openmp else 1,
                    }
                )

                text = render_nse(case, self.manifest, profile)

                self.assertIn('boundary_x_min = "reflective"', text)
                self.assertIn('boundary_x_max = "reflective"', text)
                self.assertIn('boundary_y_min = "reflective"', text)
                self.assertIn('boundary_y_max = "reflective"', text)
                self.assertIn('boundary_z_min = "reflective"', text)
                self.assertIn('boundary_z_max = "reflective"', text)
                self.assertNotIn("boundary_condition =", text)
                self.assertNotIn("boundary_x_min_reference_rho", text)

    def test_rejects_reference_state_on_reflective_face(self) -> None:
        case = self.case()
        case["numerics"].pop("boundary_condition")
        case["boundary"] = self.x_non_reflecting_boundary()
        case["boundary"]["faces"]["x_min"]["type"] = "reflective"

        with self.assertRaisesRegex(
            CaseInputError, "reference_state is only valid"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_shock_turbulence_and_dirichlet_driver(self) -> None:
        text = render_nse(
            self.shock_turbulence_case(), self.manifest, "cpu_mpi"
        )

        self.assertIn(
            'initial_condition = "shock_turbulence_interaction"', text
        )
        self.assertIn('planar_shock_direction = "positive_x"', text)
        self.assertIn("planar_shock_mach = 1.5", text)
        self.assertIn(
            "planar_shock_downstream_rho = 1.8620689655172413", text
        )
        self.assertIn(
            "planar_shock_downstream_u = 0.69444444444444442", text
        )
        self.assertIn('boundary_x_min = "dirichlet"', text)
        self.assertIn(
            "boundary_x_min_reference_rho = 1.8620689655172413", text
        )
        self.assertIn("boundary_x_min_reference_u = 0.69444444444444442", text)
        self.assertIn('boundary_x_max = "non_reflecting"', text)
        self.assertIn("boundary_x_max_reference_rho = 1", text)
        self.assertIn("imported_turbulence_background_rho = 1", text)

    def test_renders_explicit_shock_downstream_state(self) -> None:
        case = self.shock_turbulence_case()
        shock = case["flow"]["planar_shock"]
        shock.pop("mach_number")
        shock["downstream"] = {
            "density": 2.0,
            "velocity": [0.8, 0.0, 0.0],
            "pressure": 2.0,
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn("planar_shock_mach = -1", text)
        self.assertIn("planar_shock_downstream_rho = 2", text)
        self.assertIn("boundary_x_min_reference_u = 0.80000000000000004", text)

    def test_renders_negative_x_shock_with_x_max_driver(self) -> None:
        case = self.shock_turbulence_case()
        case["flow"]["planar_shock"]["propagation_direction"] = "negative_x"
        case["flow"]["planar_shock"]["position"] = case["grid"]["x_min"] + 6 * (
            (case["grid"]["x_max"] - case["grid"]["x_min"])
            / case["grid"]["nx"]
        )
        case["flow"]["imported_turbulence"]["x_start"] = case["grid"]["x_min"]
        case["boundary"]["faces"]["x_min"] = {
            "type": "non_reflecting",
            "reference_state": "shock_pre",
        }
        case["boundary"]["faces"]["x_max"] = {
            "type": "dirichlet",
            "reference_state": "shock_post",
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('planar_shock_direction = "negative_x"', text)
        self.assertIn('boundary_x_max = "dirichlet"', text)
        self.assertIn("planar_shock_downstream_u = -0.69444444444444442", text)
        self.assertIn("boundary_x_max_reference_u = -0.69444444444444442", text)

    def test_rejects_non_shock_mach_number(self) -> None:
        case = self.shock_turbulence_case()
        case["flow"]["planar_shock"]["mach_number"] = 1.0

        with self.assertRaisesRegex(CaseInputError, "must be greater than 1"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_background_different_from_upstream(self) -> None:
        case = self.shock_turbulence_case()
        case["flow"]["imported_turbulence"]["background"] = {
            "density": 2.0,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 1.0,
        }

        with self.assertRaisesRegex(CaseInputError, "must equal"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_without_dirichlet_driver(self) -> None:
        case = self.shock_turbulence_case()
        case["boundary"]["faces"]["x_min"] = {
            "type": "non_reflecting",
            "reference_state": "shock_post",
        }

        with self.assertRaisesRegex(CaseInputError, "requires.*DIRICHLET"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_driver_state_mismatch(self) -> None:
        case = self.shock_turbulence_case()
        case["boundary"]["reference_states"]["shock_post"] = {
            "density": 1.1,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 1.0,
        }

        with self.assertRaisesRegex(CaseInputError, "reference state to equal"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_reflective_shock_tube_turbulence_case(self) -> None:
        text = render_nse(
            self.shock_tube_turbulence_case(), self.manifest, "cpu_mpi"
        )

        self.assertIn(
            'initial_condition = "shock_tube_turbulence_interaction"', text
        )
        self.assertIn(
            "shock_tube_diaphragm_position = 1.5707963267948966", text
        )
        self.assertIn("shock_tube_driver_p = 3.5714285714285716", text)
        self.assertIn("shock_tube_driven_p = 0.7142857142857143", text)
        self.assertIn('boundary_x_min = "reflective"', text)
        self.assertIn('boundary_x_max = "non_reflecting"', text)
        self.assertIn("boundary_x_max_reference_rho = 1", text)
        self.assertIn("imported_turbulence_background_rho = 1", text)

    def test_rejects_shock_tube_driver_pressure_not_above_driven(self) -> None:
        case = self.shock_tube_turbulence_case()
        case["flow"]["shock_tube"]["driver"]["pressure"] = 0.5

        with self.assertRaisesRegex(CaseInputError, "pressure must be greater"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_tube_driver_normal_velocity(self) -> None:
        case = self.shock_tube_turbulence_case()
        case["flow"]["shock_tube"]["driver"]["velocity"] = [0.1, 0.0, 0.0]

        with self.assertRaisesRegex(CaseInputError, "x component must be zero"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_tube_without_reflective_closed_end(self) -> None:
        case = self.shock_tube_turbulence_case()
        case["boundary"]["faces"]["x_min"] = {
            "type": "non_reflecting",
            "reference_state": "driven",
        }

        with self.assertRaisesRegex(CaseInputError, "x_min.type=REFLECTIVE"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_tube_wrong_outflow_reference(self) -> None:
        case = self.shock_tube_turbulence_case()
        case["boundary"]["reference_states"]["wrong"] = {
            "density": 2.0,
            "velocity": [0.0, 0.0, 0.0],
            "pressure": 1.0,
        }
        case["boundary"]["faces"]["x_max"]["reference_state"] = "wrong"

        with self.assertRaisesRegex(CaseInputError, "reference state to equal"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_shock_tube_diaphragm_inside_turbulence(self) -> None:
        case = self.shock_tube_turbulence_case()
        case["flow"]["shock_tube"]["diaphragm_position"] = case["flow"][
            "imported_turbulence"
        ]["x_start"] + (
            (case["grid"]["x_max"] - case["grid"]["x_min"])
            / case["grid"]["nx"]
        )

        with self.assertRaisesRegex(CaseInputError, "must not lie inside"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_second_order_keep_selection(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "KEEP2"

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "keep2"', text)

    def test_renders_imported_turbulence_embed(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "mode": "embed",
                "x_start": 2.0,
                "blend_cells": 4,
                "velocity_offset": [0.5, 0.0, 0.0],
                "background": {
                    "density": 1.0,
                    "velocity": [0.5, 0.0, 0.0],
                    "pressure": 0.75,
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('initial_condition = "imported_turbulence"', text)
        self.assertIn(
            'imported_turbulence_file = "initial_data/turbulence.slf"', text
        )
        self.assertIn('imported_turbulence_mode = "embed"', text)
        self.assertIn("imported_turbulence_x_start = 2", text)
        self.assertIn("imported_turbulence_blend_cells = 4", text)
        self.assertIn("imported_turbulence_velocity_offset_x = 0.5", text)
        self.assertIn("imported_turbulence_background_u = 0.5", text)
        self.assertIn("imported_turbulence_background_p = 0.75", text)

    def test_resolves_imported_turbulence_file_from_case_directory(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "blend_cells": 0,
            },
        }
        runtime_root = Path.cwd() / "portable_runtime"
        case_dir = runtime_root / "cases" / "case0001"

        text = render_nse(
            case,
            self.manifest,
            "cpu_mpi",
            case_dir=case_dir,
            runtime_root=runtime_root,
        )

        self.assertIn(
            'imported_turbulence_file = '
            '"initial_data/turbulence.slf"',
            text,
        )

    def test_keeps_external_imported_turbulence_file_absolute(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "../shared/turbulence.slf",
                "blend_cells": 0,
            },
        }
        runtime_root = Path.cwd() / "portable_runtime"
        case_dir = runtime_root / "cases" / "case0001"
        external_file = (case_dir / "../shared/turbulence.slf").resolve().as_posix()

        text = render_nse(
            case,
            self.manifest,
            "cpu_mpi",
            case_dir=case_dir,
            runtime_root=runtime_root,
        )

        self.assertIn(
            f'imported_turbulence_file = "{external_file}"',
            text,
        )

    def test_rejects_tile_with_blending(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "turbulence_tile",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "blend_cells": 2,
            },
        }

        with self.assertRaisesRegex(CaseInputError, "requires blend_cells: 0"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_invalid_imported_turbulence_velocity(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "imported_turbulence",
            "imported_turbulence": {
                "file": "initial_data/turbulence.slf",
                "velocity_offset": [0.5, 0.0],
            },
        }

        with self.assertRaisesRegex(CaseInputError, "exactly three numbers"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_weno5z_roe_selection_for_cpu(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "WENO5Z_ROE"

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "weno5z_roe"', text)

    def test_renders_hybrid_keep_weno_selection(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "HYBRID"
        case["numerics"]["hybrid"] = {
            "smooth_scheme": "KEEP6",
            "shock_scheme": "WENO5Z_ROE",
            "sensor": "DUCROS_PRESSURE",
            "sensor_onset": 0.02,
            "sensor_full": 0.15,
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('convective_scheme = "hybrid"', text)
        self.assertIn('hybrid_smooth_scheme = "keep6"', text)
        self.assertIn('hybrid_shock_scheme = "weno5z_roe"', text)
        self.assertIn('hybrid_sensor = "ducros_pressure"', text)
        self.assertIn("hybrid_sensor_onset = 0.02", text)
        self.assertIn("hybrid_sensor_full = 0.14999999999999999", text)

    def test_rejects_invalid_hybrid_sensor_thresholds(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "hybrid"
        case["numerics"]["hybrid"] = {
            "sensor_onset": 0.10,
            "sensor_full": 0.05,
        }

        with self.assertRaisesRegex(CaseInputError, "sensor_onset < sensor_full"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_keep_without_order_suffix(self) -> None:
        case = self.case()
        case["numerics"]["convective_scheme"] = "KEEP"

        with self.assertRaisesRegex(CaseInputError, "KEEP2 or KEEP6"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_removed_convective_order_selector(self) -> None:
        case = self.case()
        case["numerics"]["convective_order"] = 2

        with self.assertRaisesRegex(CaseInputError, "no longer supported"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_legacy_flux_selector(self) -> None:
        case = self.case()
        case["numerics"]["flux"] = "KEEP6"

        with self.assertRaisesRegex(
            CaseInputError, r"numerics\.flux.*numerics\.convective_scheme"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_unused_reconstruction_selector(self) -> None:
        case = self.case()
        case["numerics"]["reconstruction"] = "hybrid"

        with self.assertRaisesRegex(
            CaseInputError,
            r"numerics\.reconstruction.*numerics\.convective_scheme",
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_accepts_nse_mpi_openmp_profile(self) -> None:
        _validate_solver_selection(self.case(), self.manifest, "cpu_mpi")

    def test_renders_distributed_spectral_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "spectrum": "johnsen",
                "random_seed": 24680,
                "rms_velocity": 0.1,
                "peak_wavenumber": 4.0,
                "integral_length": 1.0,
                "kolmogorov_length": 0.02,
                "dealias_fraction": 2.0 / 3.0,
                "isotropy_mode": "projected_shell",
                "isotropy_k_cutoff": 2.5,
                "isotropy_tolerance": 1.0e-8,
                "isotropy_max_iterations": 80,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('initial_condition = "hit_spectral"', text)
        self.assertIn('hit_spectrum = "johnsen"', text)
        self.assertIn("hit_seed = 24680", text)
        self.assertIn("hit_rms_velocity = 0.1", text)
        self.assertIn("hit_peak_wavenumber = 4", text)
        self.assertIn('hit_isotropy_mode = "projected_shell"', text)
        self.assertIn("hit_isotropy_k_cutoff = 2.5", text)
        self.assertIn("hit_isotropy_tolerance = 1e-08", text)
        self.assertIn("hit_isotropy_max_iterations = 80", text)
        _validate_solver_selection(
            case, self.manifest, "cpu_mpi_2decomp_fftw"
        )

    def test_rejects_hit_with_profile_without_initial_fft_backend(self) -> None:
        case = self.case()
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.1,
                "turbulent_reynolds_number": 30.0,
                "random_seed": 13579,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                },
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "initial-condition FFT backend"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_renders_target_driven_johnsen_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.3,
                "turbulent_reynolds_number": 40.0,
                "random_seed": 24680,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 2.0,
                        "length_scale_ratio": 2.5,
                    },
                    "pope": {"integral_length": 1.0},
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")
        derived = derive_nse_hit_transport(0.3, 40.0, 2.0)

        self.assertIn('hit_spectrum = "johnsen"', text)
        self.assertIn("hit_turbulent_mach = 0.29999999999999999", text)
        self.assertIn("hit_turbulent_reynolds = 40", text)
        self.assertIn("hit_johnsen_length_scale_ratio = 2.5", text)
        self.assertIn("hit_peak_wavenumber = 2.5", text)
        self.assertIn(
            f"reynolds = {format(derived['solver_reynolds'], '.17g')}", text
        )

    def test_renders_target_driven_pope_hit_input(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "taylor_reynolds_number": 30.0,
                "spectrum": {
                    "type": "pope",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                    "pope": {
                        "integral_length": 1.0,
                        "energy_constant": 1.5,
                        "large_scale_constant": 6.78,
                        "dissipation_constant": 0.4,
                        "large_scale_exponent": 2.0,
                        "dissipation_exponent": 5.2,
                    },
                },
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")
        derived = derive_nse_hit_transport(0.5, 30.0, 1.0)

        self.assertIn('hit_spectrum = "pope"', text)
        self.assertIn("hit_pope_energy_constant = 1.5", text)
        self.assertIn("hit_pope_large_scale_constant = 6.7800000000000002", text)
        self.assertIn("hit_pope_dissipation_exponent = 5.2000000000000002", text)
        self.assertIn(
            "hit_kolmogorov_length = "
            f"{format(derived['kolmogorov_length'], '.17g')}",
            text,
        )

    def test_rejects_unknown_hit_spectrum_type(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "turbulent_reynolds_number": 30.0,
                "spectrum": {"type": "poppe"},
            },
        }

        with self.assertRaisesRegex(CaseInputError, "spectrum.type"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_hit_target_missing_reynolds_number(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["flow"] = {
            "type": "hit",
            "hit": {
                "turbulent_mach_number": 0.5,
                "spectrum": {
                    "type": "johnsen",
                    "johnsen": {
                        "characteristic_length": 1.0,
                        "length_scale_ratio": 2.0,
                    },
                },
            },
        }

        with self.assertRaisesRegex(CaseInputError, "turbulent_reynolds_number"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_renders_petersen_livescu_forcing(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "low_wavenumber",
                "fft_backend": "2decomp_fftw",
                "k_cutoff": 2.5,
                "target_dissipation": 0.1,
                "dilatational_ratio": 0.25,
                "denominator_floor": 1.0e-14,
                "max_coefficient": 20.0,
                "report_interval": 50,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)
        self.assertIn('forcing_spectrum = "low_wavenumber"', text)
        self.assertIn('forcing_fft_backend = "2decomp_fftw"', text)
        self.assertIn("forcing_k_cutoff = 2.5", text)
        self.assertIn("forcing_target_dissipation = 0.1", text)
        self.assertIn("forcing_dilatational_ratio = 0.25", text)
        self.assertIn("forcing_report_interval = 50", text)

    def test_rejects_forcing_with_non_periodic_boundary(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["numerics"].pop("boundary_condition")
        case["boundary"] = self.x_non_reflecting_boundary()
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {"target_dissipation": 0.1},
        }

        with self.assertRaisesRegex(CaseInputError, "periodic boundaries"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_normalizes_human_readable_forcing_selectors(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "Petersen-Livescu",
            "petersen_livescu": {
                "spectrum": "Low Wavenumber",
                "fft_backend": "2decomp-fftw",
                "target_dissipation": 0.1,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)
        self.assertIn('forcing_spectrum = "low_wavenumber"', text)
        self.assertIn('forcing_fft_backend = "2decomp_fftw"', text)

    def test_normalizes_full_wavenumber_forcing_alias(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "full_wavenumber",
                "target_dissipation": 0.1,
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_spectrum = "full_spectrum"', text)

    def test_rejects_unknown_forcing_spectrum(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "spectrum": "full_wavenumbar",
                "target_dissipation": 0.1,
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "forcing.petersen_livescu.spectrum"
        ):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_forcing_with_profile_without_fft_backend(self) -> None:
        case = self.case()
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {"target_dissipation": 0.1},
        }

        with self.assertRaisesRegex(
            CaseInputError, "compatible staged profile"
        ):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_accepts_legacy_forcing_scheme_alias(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "scheme": "petersen_livescu",
            "target_dissipation": 0.1,
        }

        text = render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

        self.assertIn('forcing_scheme = "petersen_livescu"', text)

    def test_ignores_inactive_petersen_livescu_settings(self) -> None:
        case = self.case()
        case["forcing"] = {
            "type": "none",
            "petersen_livescu": {
                "target_dissipation": 0.1,
                "fft_backend": "2decomp_fftw",
            },
        }

        text = render_nse(case, self.manifest, "cpu_mpi")

        self.assertIn('forcing_scheme = "none"', text)
        self.assertNotIn("forcing_target_dissipation", text)

    def test_rejects_unknown_forcing_type(self) -> None:
        case = self.case()
        case["forcing"] = {"type": "petersen_livecu"}

        with self.assertRaisesRegex(CaseInputError, "unknown forcing.type"):
            render_nse(case, self.manifest, "cpu_mpi")

    def test_rejects_misspelled_type_specific_forcing_key(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "petersen_livescu": {
                "target_disipation": 0.1,
            },
        }

        with self.assertRaisesRegex(
            CaseInputError, "unknown forcing.petersen_livescu key"
        ):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_rejects_mixed_nested_and_legacy_forcing_settings(self) -> None:
        case = self.case()
        case["solver"]["profile"] = "cpu_mpi_2decomp_fftw"
        case["forcing"] = {
            "type": "petersen_livescu",
            "target_dissipation": 0.1,
            "petersen_livescu": {"target_dissipation": 0.2},
        }

        with self.assertRaisesRegex(CaseInputError, "mixes the new"):
            render_nse(case, self.manifest, "cpu_mpi_2decomp_fftw")

    def test_renders_nse_single_gpu_input(self) -> None:
        case = self.case(mpi_processes=1)
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
                "cuda_device": 2,
            }
        )
        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('backend = "cuda"', text)
        self.assertIn("use_mpi = .false.", text)
        self.assertIn("use_openmp = .false.", text)
        self.assertIn("cuda_device = 2", text)
        _validate_solver_selection(case, self.manifest, "cuda_single")

    def test_renders_weno5z_roe_for_cuda_profile(self) -> None:
        case = self.case(mpi_processes=1)
        case["solver"].update(
            {
                "profile": "cuda_single",
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )
        case["numerics"]["convective_scheme"] = "WENO5Z_ROE"

        text = render_nse(case, self.manifest, "cuda_single")

        self.assertIn('convective_scheme = "weno5z_roe"', text)

    def test_renders_nse_multi_gpu_input(self) -> None:
        case = self.case(mpi_processes=2)
        case["solver"].update(
            {
                "profile": "cuda_mpi",
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": True,
                "omp_threads": 1,
            }
        )

        text = render_nse(case, self.manifest, "cuda_mpi")

        self.assertIn('backend = "cuda_mpi"', text)
        self.assertIn("use_mpi = .true.", text)
        self.assertIn("use_openmp = .false.", text)
        _validate_solver_selection(case, self.manifest, "cuda_mpi")

    def test_rejects_too_few_nse_mpi_processes(self) -> None:
        with self.assertRaisesRegex(CaseInputError, "at least 2"):
            _validate_solver_selection(
                self.case(mpi_processes=1), self.manifest, "cpu_mpi"
            )


if __name__ == "__main__":
    unittest.main()
