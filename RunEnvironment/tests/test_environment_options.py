from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from environment_options import (  # noqa: E402
    OptionCatalogError,
    load_option_catalogs,
    resolve_design_options,
)
from yaml_support import load_yaml  # noqa: E402


BUILTIN = SCRIPT_DIR / "environment_options.yaml"


class EnvironmentOptionTests(unittest.TestCase):
    def test_supplied_design_templates_resolve(self) -> None:
        expected = {
            "environment.yaml": ("cpu_mpi_fftw", "workstation", False),
            "environment.hpc.yaml": ("cpu_mpi_fftw", "hpc_slurm", True),
        }
        for name, values in expected.items():
            with self.subTest(name=name):
                design_path = SCRIPT_DIR / name
                resolved, _ = resolve_design_options(
                    load_yaml(design_path), design_path, BUILTIN
                )
                self.assertEqual(resolved["model"]["profile"], values[0])
                self.assertEqual(resolved["target"]["type"], values[1])
                self.assertEqual(resolved["scheduler"]["enabled"], values[2])

    def test_builtin_selection_and_inline_override(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "source": "mozart_nas",
                "model": "gpe_cpu_mpi_fftw",
                "case": "gpe_quantum_taylor_green",
                "execution": "mpi4_release",
            },
            "execution": {"processes": 6},
        }
        resolved, state = resolve_design_options(
            design, SCRIPT_DIR / "environment.yaml", BUILTIN
        )

        self.assertEqual(
            resolved["source"]["framework_root"],
            r"\\Mozart\share\FrameWork",
        )
        self.assertEqual(resolved["model"]["profile"], "cpu_mpi_fftw")
        self.assertEqual(resolved["execution"]["processes"], 6)
        self.assertEqual(resolved["execution"]["configuration"], "Release")
        self.assertEqual(state["selections"]["model"]["id"], "gpe_cpu_mpi_fftw")

    def test_legacy_inline_design_needs_no_selection(self) -> None:
        design = {
            "schema_version": 1,
            "source": {"framework_root": "../.."},
            "model": {"name": "nse", "profile": "cpu_mpi"},
        }
        resolved, state = resolve_design_options(
            design, SCRIPT_DIR / "environment.yaml", BUILTIN
        )

        self.assertEqual(resolved, design)
        self.assertEqual(state["selections"], {})

    def test_custom_catalog_adds_an_option(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog_path = root / "local.yaml"
            catalog_path.write_text(
                "\n".join(
                    [
                        "schema_version: 1",
                        "catalog_id: local_test",
                        "options:",
                        "  execution:",
                        "    mpi16_release:",
                        "      label: MPI 16",
                        "      values:",
                        "        configuration: Release",
                        "        processes: 16",
                        "        omp_threads: 1",
                        "        parallel_jobs: 16",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            design_path = root / "environment.yaml"
            design_path.write_text(
                "\n".join(
                    [
                        "schema_version: 1",
                        "option_catalogs:",
                        "  - local.yaml",
                        "select:",
                        "  execution: mpi16_release",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            design = load_yaml(design_path)
            resolved, state = resolve_design_options(
                design, design_path, BUILTIN
            )

            self.assertEqual(resolved["execution"]["processes"], 16)
            self.assertIn("mpi16_release", state["options"]["execution"])

    def test_duplicate_option_requires_explicit_override(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog_path = root / "duplicate.yaml"
            catalog_path.write_text(
                "\n".join(
                    [
                        "schema_version: 1",
                        "catalog_id: duplicate_test",
                        "options:",
                        "  archive:",
                        "    none:",
                        "      values:",
                        "        enabled: true",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            design = {
                "schema_version": 1,
                "option_catalogs": [str(catalog_path)],
            }
            with self.assertRaises(OptionCatalogError):
                load_option_catalogs(design, root / "environment.yaml", BUILTIN)

    def test_incompatible_model_and_case_are_rejected(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "model": "nse_cpu_mpi",
                "case": "gpe_quantum_taylor_green",
            },
        }
        with self.assertRaises(OptionCatalogError):
            resolve_design_options(
                design, SCRIPT_DIR / "environment.yaml", BUILTIN
            )

    def test_single_gpu_rejects_multiple_processes(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "model": "gpe_cuda_single",
                "execution": "mpi8_release",
            },
        }
        with self.assertRaises(OptionCatalogError):
            resolve_design_options(
                design, SCRIPT_DIR / "environment.yaml", BUILTIN
            )

    def test_single_gpu_accepts_serial_execution(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "model": "gpe_cuda_single",
                "execution": "serial_release",
            },
        }
        resolved, _ = resolve_design_options(
            design, SCRIPT_DIR / "environment.yaml", BUILTIN
        )

        self.assertEqual(resolved["model"]["profile"], "cuda_single")
        self.assertEqual(resolved["execution"]["processes"], 1)


if __name__ == "__main__":
    unittest.main()
