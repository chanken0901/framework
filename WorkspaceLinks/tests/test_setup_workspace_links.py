from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from setup_workspace_links import (  # noqa: E402
    DesignIntegration,
    LinkConfig,
    LinkSpec,
    WorkspaceLinkError,
    _apply,
    _design_status,
    _inspect_links,
    _remove,
)


class WorkspaceLinkTests(unittest.TestCase):
    def _config(self, root: Path) -> LinkConfig:
        framework = root / "library" / "FrameWork"
        framework.mkdir(parents=True)
        return LinkConfig(
            source=root / "workspace_links.yaml",
            execution_root=root / "ResearchRuns",
            link_type="auto",
            designs=DesignIntegration(
                legacy_path=root / "ResearchDesigns",
                integrated_path=root / "ResearchRuns" / "Designs",
                keep_legacy_link=True,
            ),
            links=(LinkSpec("framework", "FrameWork", framework),),
        )

    def test_apply_moves_designs_and_remove_only_unlinks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self._config(root)
            config.designs.legacy_path.mkdir()
            design_file = config.designs.legacy_path / "gpe_case0001.yaml"
            design_file.write_text("schema_version: 1\n", encoding="utf-8")

            _apply(config, replace_links=False, dry_run=False)

            integrated_file = config.designs.integrated_path / design_file.name
            self.assertTrue(integrated_file.is_file())
            self.assertEqual(_design_status(config), "ready")
            self.assertTrue(
                all(state.status == "ok" for state in _inspect_links(config))
            )

            _remove(config, dry_run=False)

            self.assertFalse(os.path.lexists(config.designs.legacy_path))
            self.assertFalse(
                os.path.lexists(config.execution_root / config.links[0].name)
            )
            self.assertTrue(integrated_file.is_file())
            self.assertTrue(config.links[0].target.is_dir())

    def test_apply_merges_non_conflicting_design_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self._config(root)
            config.designs.legacy_path.mkdir()
            config.designs.integrated_path.mkdir(parents=True)
            (config.designs.legacy_path / "gpe.yaml").write_text("", encoding="utf-8")
            (config.designs.integrated_path / "nse.yaml").write_text("", encoding="utf-8")

            _apply(config, replace_links=False, dry_run=False)

            self.assertTrue((config.designs.integrated_path / "gpe.yaml").is_file())
            self.assertTrue((config.designs.integrated_path / "nse.yaml").is_file())
            self.assertEqual(_design_status(config), "ready")

    def test_apply_refuses_conflicting_design_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self._config(root)
            config.designs.legacy_path.mkdir()
            config.designs.integrated_path.mkdir(parents=True)
            (config.designs.legacy_path / "case.yaml").write_text("old", encoding="utf-8")
            (config.designs.integrated_path / "case.yaml").write_text(
                "new", encoding="utf-8"
            )

            with self.assertRaises(WorkspaceLinkError):
                _apply(config, replace_links=False, dry_run=False)

            self.assertEqual(
                (config.designs.legacy_path / "case.yaml").read_text(encoding="utf-8"),
                "old",
            )
            self.assertEqual(
                (config.designs.integrated_path / "case.yaml").read_text(
                    encoding="utf-8"
                ),
                "new",
            )

    def test_apply_refuses_to_replace_a_real_framework_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = self._config(root)
            config.execution_root.mkdir()
            (config.execution_root / "FrameWork").mkdir()

            with self.assertRaises(WorkspaceLinkError):
                _apply(config, replace_links=True, dry_run=False)


if __name__ == "__main__":
    unittest.main()
