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
            "environment.gpe.yaml": ("gpe", "workstation", False),
            "environment.nse.yaml": ("nse", "workstation", False),
            "environment.hpc.yaml": ("gpe", "hpc_slurm", True),
        }
        for name, values in expected.items():
            with self.subTest(name=name):
                design_path = SCRIPT_DIR / name
                resolved, _ = resolve_design_options(
                    load_yaml(design_path), design_path, BUILTIN
                )
                self.assertEqual(resolved["model"]["name"], values[0])
                self.assertEqual(resolved["target"]["type"], values[1])
                self.assertEqual(resolved["scheduler"]["enabled"], values[2])

    def test_builtin_selection_and_inline_override(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "source": "mozart_nas",
                "model": "gpe",
                "case": "gpe_quantum_taylor_green",
                "execution": "release",
            },
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
            },
        }
        resolved, state = resolve_design_options(
            design, SCRIPT_DIR / "environment.gpe.yaml", BUILTIN
        )

        self.assertEqual(
            resolved["source"]["framework_root"],
            r"\\Mozart\share\FrameWork",
        )
        self.assertEqual(resolved["model"]["name"], "gpe")
        self.assertEqual(resolved["execution"]["configuration"], "Release")
        self.assertTrue(resolved["parallel"]["use_mpi"])
        self.assertEqual(state["selections"]["model"]["id"], "gpe")

    def test_legacy_composite_model_preserves_profile_override(self) -> None:
        design = {
            "schema_version": 1,
            "select": {"model": "gpe_cpu_mpi_dft"},
            "parallel": {
                "use_mpi": True,
                "use_openmp": False,
                "use_cuda": False,
            },
        }

        resolved, state = resolve_design_options(
            design, SCRIPT_DIR / "environment.gpe.yaml", BUILTIN
        )

        self.assertEqual(resolved["model"]["name"], "gpe")
        self.assertEqual(resolved["solver"]["profile"], "cpu_mpi_dft")
        self.assertEqual(state["selections"]["model"]["id"], "gpe")

    def test_legacy_inline_design_needs_no_selection(self) -> None:
        design = {
            "schema_version": 1,
            "source": {"framework_root": "../.."},
            "model": {"name": "nse", "profile": "cpu_mpi"},
        }
        resolved, state = resolve_design_options(
            design, SCRIPT_DIR / "environment.gpe.yaml", BUILTIN
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
                        "    release_large_build:",
                        "      label: Release large build",
                        "      values:",
                        "        configuration: Release",
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
                        "  execution: release_large_build",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            design = load_yaml(design_path)
            resolved, state = resolve_design_options(
                design, design_path, BUILTIN
            )

            self.assertEqual(resolved["execution"]["parallel_jobs"], 16)
            self.assertIn("release_large_build", state["options"]["execution"])

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
                "model": "nse",
                "case": "gpe_quantum_taylor_green",
            },
        }
        with self.assertRaises(OptionCatalogError):
            resolve_design_options(
                design, SCRIPT_DIR / "environment.gpe.yaml", BUILTIN
            )

    def test_single_gpu_uses_release_execution_without_parallel_count(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "model": "gpe",
                "execution": "release",
            },
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
            },
        }
        resolved, _ = resolve_design_options(
            design, SCRIPT_DIR / "environment.gpe.yaml", BUILTIN
        )

        self.assertEqual(resolved["model"]["name"], "gpe")
        self.assertNotIn("processes", resolved["execution"])
        self.assertTrue(resolved["parallel"]["use_cuda"])

    def test_nse_single_gpu_accepts_cuda_parallel_flags(self) -> None:
        design = {
            "schema_version": 1,
            "select": {
                "model": "nse",
                "case": "nse_case",
                "execution": "release",
            },
            "parallel": {
                "use_mpi": False,
                "use_openmp": False,
                "use_cuda": True,
            },
        }
        resolved, _ = resolve_design_options(
            design, SCRIPT_DIR / "environment.nse.yaml", BUILTIN
        )

        self.assertEqual(resolved["model"]["name"], "nse")
        self.assertNotIn("processes", resolved["execution"])

    def test_stage_four_multicomponent_environment_resolves(self) -> None:
        design_path = (
            SCRIPT_DIR / "environment.nse_multicomponent.viscous.yaml"
        )
        resolved, _ = resolve_design_options(
            load_yaml(design_path), design_path, BUILTIN
        )

        self.assertEqual(resolved["model"]["name"], "nse_multicomponent")
        self.assertEqual(resolved["solver"]["profile"], "cpu_serial_viscous")
        self.assertEqual(
            resolved["case"]["template"],
            "ScriptLibrary/RunEnvironment/case_templates/"
            "nse_multicomponent_viscous.yaml",
        )
        self.assertEqual(
            set(resolved["case"]["extension_templates"]),
            {"multicomponent", "thermodynamics", "transport"},
        )
        self.assertFalse(resolved["parallel"]["use_mpi"])
        self.assertFalse(resolved["parallel"]["use_cuda"])

    def test_stage_five_multicomponent_environment_resolves(self) -> None:
        design_path = (
            SCRIPT_DIR / "environment.nse_multicomponent.reactor.yaml"
        )
        resolved, _ = resolve_design_options(
            load_yaml(design_path), design_path, BUILTIN
        )

        self.assertEqual(resolved["model"]["name"], "nse_multicomponent")
        self.assertEqual(resolved["solver"]["profile"], "cpu_serial_reactor")
        self.assertEqual(
            set(resolved["case"]["extension_templates"]),
            {"multicomponent", "thermodynamics", "chemistry"},
        )
        self.assertFalse(resolved["parallel"]["use_mpi"])
        self.assertFalse(resolved["parallel"]["use_cuda"])


if __name__ == "__main__":
    unittest.main()
