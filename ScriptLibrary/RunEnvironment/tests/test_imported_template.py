"""Regression: generated templates expose mode-specific required inputs."""
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from prepare_environment import prepare
from case_input import CaseInputError, render_nse
from yaml_support import load_yaml


class ImportedTemplateTests(unittest.TestCase):
    def test_generated_template_can_switch_modes(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = Namespace(design=str(SCRIPT_DIR / "environment.nse.yaml"),
                framework_root=str(SCRIPT_DIR.parents[1]), output=str(Path(temporary) / "case"),
                model=None, profile=None, case_id=None, overwrite=False, archive=False,
                archive_format=None, dry_run=False)
            with patch("prepare_environment.sync_environment_case"):
                generated = prepare(args)
            case = load_yaml(generated / "cases/case0001/case.yaml")
            manifest = load_yaml(generated / "SolverLibrary/NSE/solver_manifest.yaml")
            imported = case["flow"]["imported_turbulence"]
            self.assertIn("x_length", imported)
            self.assertIsNone(imported["x_length"])
            case["flow"]["type"] = "imported_turbulence"
            imported["blend_cells"] = 0
            for mode in ("embed", "tile"):
                imported["mode"] = mode
                self.assertNotIn("imported_turbulence_x_length", render_nse(case, manifest, "cpu_mpi"))
            imported["mode"] = "periodic_embed"
            with self.assertRaisesRegex(CaseInputError, "x_length is required"):
                render_nse(case, manifest, "cpu_mpi")
            imported["x_length"] = 4*(case["grid"]["x_max"]-case["grid"]["x_min"])/case["grid"]["nx"]
            self.assertIn("imported_turbulence_x_length", render_nse(case, manifest, "cpu_mpi"))
