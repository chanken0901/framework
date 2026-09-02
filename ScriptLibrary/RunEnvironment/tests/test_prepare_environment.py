from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from prepare_environment import (  # noqa: E402
    EnvironmentError,
    _active_fortran_source,
    _compatible_profile_names,
    _global_case_index_path,
    _inspect_dependencies,
    _numbered_case_destination,
    _parallel_features,
    _profile_preprocessor_defines,
    _profile_is_explicit,
    _resolve_solver_profile,
    _requested_case_destination,
    _selected_solver_files,
    _slurm_script,
)
from yaml_support import load_yaml  # noqa: E402


class NumberedCaseDestinationTests(unittest.TestCase):
    def test_uses_first_available_model_case_number(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "gpe_case0001").mkdir()
            (root / "gpe_case0002.tar.gz").touch()

            output, case_id = _numbered_case_destination(root, "gpe")

            self.assertEqual(output, root / "gpe_case0003")
            self.assertEqual(case_id, "case0003")

    def test_supports_an_independent_model_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "gpe_case0001").mkdir()

            output, case_id = _numbered_case_destination(root, "nse")

            self.assertEqual(output, root / "nse_case0001")
            self.assertEqual(case_id, "case0001")

    def test_selects_an_existing_case_number_for_overwrite(self) -> None:
        root = Path("C:/ResearchRuns")

        output, case_id = _requested_case_destination(root, "nse", "case0015")

        self.assertEqual(output, root / "nse_case0015")
        self.assertEqual(case_id, "case0015")

    def test_normalizes_a_numeric_overwrite_case_id(self) -> None:
        root = Path("C:/ResearchRuns")

        output, case_id = _requested_case_destination(root, "gpe", "15")

        self.assertEqual(output, root / "gpe_case0015")
        self.assertEqual(case_id, "case0015")

    def test_rejects_an_invalid_overwrite_case_id(self) -> None:
        with self.assertRaisesRegex(EnvironmentError, "--case-id"):
            _requested_case_destination(Path("C:/ResearchRuns"), "nse", "latest")

    def test_global_index_must_be_outside_case_environment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "gpe_case0001"

            with self.assertRaises(EnvironmentError):
                _global_case_index_path(
                    {"case_index": str(output / "cases.csv")},
                    root / "environment.yaml",
                    output,
                )


class SlurmScriptTests(unittest.TestCase):
    def test_generated_job_runs_existing_executable_only(self) -> None:
        script = _slurm_script(
            {
                "job_name": "gpe-case",
                "nodes": 1,
                "tasks_per_node": 4,
                "modules": ["gcc", "openmpi"],
            },
            include_tests=False,
            use_mpi=True,
            use_openmp=True,
        )

        self.assertIn("python3 tools/run_case.py --run", script)
        self.assertIn('--processes "${SLURM_NTASKS}"', script)
        self.assertIn('--omp-threads "${SLURM_CPUS_PER_TASK}"', script)
        self.assertNotIn("#SBATCH --ntasks-per-node", script)
        self.assertNotIn("#SBATCH --cpus-per-task", script)
        self.assertNotIn("--build", script)
        self.assertNotIn("--all", script)


class ParallelFeatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        framework_root = SCRIPT_DIR.parents[1]
        cls.gpe_manifest = load_yaml(
            framework_root
            / "SolverLibrary"
            / "GPE"
            / "gp3d"
            / "solver_manifest.yaml"
        )
        cls.nse_manifest = load_yaml(
            framework_root / "SolverLibrary" / "NSE" / "solver_manifest.yaml"
        )

    def test_accepts_mpi_openmp_capability_without_counts(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
            }
        }

        self.assertEqual(
            _parallel_features(design, self.gpe_manifest, "cpu_mpi_fftw"),
            (True, True, False),
        )

    def test_rejects_cuda_flag_for_cpu_profile(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": True,
            }
        }

        with self.assertRaisesRegex(EnvironmentError, "use_cuda"):
            _parallel_features(design, self.gpe_manifest, "cpu_mpi_fftw")

    def test_gpe_profile_is_derived_from_mpi_and_cuda(self) -> None:
        expected = {
            (False, False): "cpu_serial_fftw",
            (True, False): "cpu_mpi_fftw",
            (False, True): "cuda_single",
            (True, True): "cuda_mpi_cufftmp",
        }
        for (use_mpi, use_cuda), profile in expected.items():
            with self.subTest(use_mpi=use_mpi, use_cuda=use_cuda):
                design = {
                    "parallel": {
                        "use_mpi": use_mpi,
                        "use_openmp": False,
                        "use_cuda": use_cuda,
                    }
                }
                resolved = _resolve_solver_profile(design, self.gpe_manifest)
                self.assertEqual(resolved[0], profile)

    def test_gpe_pencil_decomposition_selects_2decomp_profile(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
                "fft_decomposition": "pencil",
            }
        }

        resolved = _resolve_solver_profile(design, self.gpe_manifest)

        self.assertEqual(resolved[0], "cpu_mpi_pencil_fftw")

    def test_gpe_pencil_decomposition_selects_cufftmp_profile(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": True,
                "fft_decomposition": "pencil",
            }
        }

        resolved = _resolve_solver_profile(design, self.gpe_manifest)

        self.assertEqual(resolved[0], "cuda_mpi_cufftmp_pencil")

    def test_gpe_pencil_decomposition_rejects_single_gpu(self) -> None:
        design = {
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
                "fft_decomposition": "pencil",
            }
        }

        with self.assertRaisesRegex(EnvironmentError, "multi-GPU"):
            _resolve_solver_profile(design, self.gpe_manifest)

    def test_explicit_slab_profile_rejects_pencil_selection(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
                "fft_decomposition": "pencil",
            },
            "solver": {"profile": "cpu_mpi_fftw"},
        }

        with self.assertRaisesRegex(EnvironmentError, "fft_decomposition"):
            _resolve_solver_profile(design, self.gpe_manifest)

    def test_nse_profile_is_derived_from_mpi_and_cuda(self) -> None:
        mpi = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": True,
                "use_cuda": False,
            }
        }
        cuda = {
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
            }
        }

        self.assertEqual(
            _resolve_solver_profile(mpi, self.nse_manifest)[0], "cpu_mpi"
        )
        self.assertEqual(
            _resolve_solver_profile(cuda, self.nse_manifest)[0], "cuda_single"
        )

    def test_automatic_nse_mpi_environment_stages_fft_profile(self) -> None:
        self.assertTrue(self.nse_manifest["runtime_profile_selection"])
        self.assertEqual(
            _compatible_profile_names(
                self.nse_manifest, "cpu_mpi", require_openmp=True
            ),
            ["cpu_mpi", "cpu_mpi_2decomp_fftw"],
        )

    def test_solver_profile_override_is_recognized_as_explicit(self) -> None:
        self.assertFalse(_profile_is_explicit({"model": {"name": "nse"}}))
        self.assertTrue(
            _profile_is_explicit(
                {
                    "model": {"name": "nse"},
                    "solver": {"profile": "cpu_mpi"},
                }
            )
        )

    def test_explicit_specialized_profile_must_match_parallel_flags(self) -> None:
        design = {
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
            },
            "solver": {"profile": "cpu_mpi_dft"},
        }

        self.assertEqual(
            _resolve_solver_profile(design, self.gpe_manifest)[0],
            "cpu_mpi_dft",
        )

    def test_explicit_profile_rejects_parallel_mismatch(self) -> None:
        design = {
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
            },
            "solver": {"profile": "cpu_mpi_dft"},
        }

        with self.assertRaisesRegex(EnvironmentError, "use_mpi"):
            _resolve_solver_profile(design, self.gpe_manifest)

    def test_nse_rejects_unsupported_cpu_serial_mode(self) -> None:
        design = {
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": False,
            }
        }

        with self.assertRaisesRegex(EnvironmentError, "no default profile"):
            _resolve_solver_profile(design, self.nse_manifest)

    def test_legacy_model_profile_is_still_accepted(self) -> None:
        design = {
            "model": {"name": "gpe", "profile": "cpu_mpi_dft"},
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
            },
        }

        self.assertEqual(
            _resolve_solver_profile(design, self.gpe_manifest)[0],
            "cpu_mpi_dft",
        )


class FortranDependencyInspectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.framework_root = SCRIPT_DIR.parents[1]
        cls.solver_root = cls.framework_root / "SolverLibrary" / "NSE"
        cls.manifest = load_yaml(cls.solver_root / "solver_manifest.yaml")

    def test_preprocessor_selects_only_the_active_use_branch(self) -> None:
        source = """\
#ifdef USE_GPU
  use mod_gpu
#else
  use mod_cpu
#endif
#if defined(USE_EXTRA)
  use mod_extra
#endif
"""

        cpu = _active_fortran_source(source, set(), "conditional.f90")
        gpu = _active_fortran_source(
            source, {"USE_GPU", "USE_EXTRA"}, "conditional.f90"
        )

        self.assertIn("use mod_cpu", cpu)
        self.assertNotIn("use mod_gpu", cpu)
        self.assertNotIn("use mod_extra", cpu)
        self.assertIn("use mod_gpu", gpu)
        self.assertIn("use mod_extra", gpu)
        self.assertNotIn("use mod_cpu", gpu)

    def test_real_nse_fft_profiles_pass_dependency_inspection(self) -> None:
        for profile_name in ("cpu_mpi_2decomp_fftw", "cuda_mpi_cufftmp"):
            with self.subTest(profile=profile_name):
                files, _ = _selected_solver_files(
                    self.solver_root, self.manifest, profile_name, False
                )
                dependencies = _inspect_dependencies(
                    self.solver_root,
                    self.manifest,
                    files,
                    _profile_preprocessor_defines(self.manifest, profile_name),
                )
                self.assertIn(
                    "src/init/mod_init_hit_spectral_2decomp.f90", dependencies
                )

    def test_active_cufftmp_use_still_requires_its_provider(self) -> None:
        files, _ = _selected_solver_files(
            self.solver_root,
            self.manifest,
            "cpu_mpi_2decomp_fftw",
            False,
        )

        with self.assertRaisesRegex(EnvironmentError, "mod_nse_cufftmp_fft"):
            _inspect_dependencies(
                self.solver_root,
                self.manifest,
                files,
                {"NSE_INIT_CUFFTMP"},
            )


class NseCaseTemplateTests(unittest.TestCase):
    def test_parallel_features_and_runtime_counts_are_template_fields(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("use_openmp: {{use_openmp}}", template)
        self.assertIn("mpi_processes: {{mpi_processes}}", template)
        self.assertIn("omp_threads: {{omp_threads}}", template)
        self.assertNotIn("profile: {{solver_profile}}", template)
        self.assertNotIn("use_mpi: {{use_mpi}}", template)
        self.assertNotIn("use_cuda: {{use_cuda}}", template)
        self.assertNotIn("use_openmp: true", template)

    def test_gpe_template_omits_profile_derived_fields(self) -> None:
        template = (
            SCRIPT_DIR
            / "case_templates"
            / "gpe_quantum_taylor_green.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("use_openmp: {{use_openmp}}", template)
        self.assertNotIn("profile: {{solver_profile}}", template)
        self.assertNotIn("use_mpi: {{use_mpi}}", template)
        self.assertNotIn("use_cuda: {{use_cuda}}", template)
        self.assertIn("parallel.fft_decomposition", template)

    def test_flow_conditions_share_one_nse_template(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("type: taylor_green", template)
        self.assertIn("taylor_green:", template)
        self.assertIn("hit:", template)

    def test_forcing_conditions_use_type_specific_mapping(self) -> None:
        template = (
            SCRIPT_DIR / "case_templates" / "nse.yaml"
        ).read_text(encoding="utf-8")

        self.assertIn("forcing:\n  # choices: none, petersen_livescu", template)
        self.assertIn("  type: none\n  petersen_livescu:\n", template)
        self.assertIn(
            "# choices: full_spectrum, low_wavenumber", template
        )
        self.assertIn("# choices: auto, 2decomp_fftw, cufft", template)
        self.assertIn("    target_dissipation: 0.1", template)
        self.assertNotIn("\n  target_dissipation: 0.1", template)


if __name__ == "__main__":
    unittest.main()
