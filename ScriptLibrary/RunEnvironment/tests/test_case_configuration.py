from __future__ import annotations

import sys
import tempfile
import unittest
from datetime import date, datetime, timedelta, timezone
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from case_configuration import (  # noqa: E402
    CaseConfigurationError,
    ResolvedCaseConfiguration,
    resolve_case_configuration,
    write_resolved_case,
)
from yaml_support import load_yaml  # noqa: E402


class CaseConfigurationTests(unittest.TestCase):
    def test_resolved_case_preserves_dates_and_nested_timestamps(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "case.yaml"
            stamp = datetime(2026, 9, 7, 12, 34, 56, 123456,
                             tzinfo=timezone(timedelta(hours=9)))
            document = {"created": stamp, "metadata": {
                "dates": [date(2026, 9, 7), datetime(2026, 9, 7, 1, 2, 3)],
                "label": "衝撃波", "enabled": True, "steps": 1000,
            }}
            resolved = ResolvedCaseConfiguration(document, (source,), {})
            output = write_resolved_case(resolved)
            reloaded = load_yaml(output)
            self.assertEqual(reloaded["created"], stamp.isoformat())
            self.assertEqual(reloaded["metadata"]["dates"],
                             ["2026-09-07", "2026-09-07T01:02:03"])
            self.assertEqual(reloaded["metadata"]["steps"], 1000)
            self.assertIs(reloaded["metadata"]["enabled"], True)
            self.assertIs(document["created"], stamp)

    def test_unquoted_yaml_timestamp_can_be_resolved_and_written(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "case.yaml"
            source.write_text("created: 2026-09-07T12:34:56+09:00\n"
                              "date: 2026-09-07\nphysics:\n  model: nse\n",
                              encoding="utf-8")
            output = write_resolved_case(resolve_case_configuration(source))
            self.assertEqual(load_yaml(output)["created"],
                             "2026-09-07T12:34:56+09:00")
            self.assertEqual(load_yaml(output)["date"], "2026-09-07")

    def test_unsupported_object_does_not_replace_existing_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "case.yaml"
            destination = Path(temporary) / "resolved_case.yaml"
            destination.write_text('{"existing": true}\n', encoding="utf-8")
            resolved = ResolvedCaseConfiguration({"bad": object()}, (source,), {})
            with self.assertRaises(TypeError):
                write_resolved_case(resolved)
            self.assertEqual(load_yaml(destination), {"existing": True})
            self.assertEqual(list(Path(temporary).glob(".*.tmp")), [])

    @staticmethod
    def _write_modular_case(root: Path, extensions: str) -> Path:
        case_path = root / "case.yaml"
        case_path.write_text(
            "\n".join(
                [
                    "schema_version: 2",
                    "case_id: case0001",
                    "physics:",
                    "  model: nse_multicomponent",
                    "extensions:",
                    extensions,
                    "",
                ]
            ),
            encoding="utf-8",
        )
        return case_path

    def test_legacy_single_file_case_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            case_path = Path(temporary) / "case.yaml"
            case_path.write_text(
                "physics:\n  model: nse\nthermodynamics:\n  gamma: 1.4\n",
                encoding="utf-8",
            )

            resolved = resolve_case_configuration(case_path)

            self.assertEqual(resolved.document["physics"]["model"], "nse")
            self.assertEqual(resolved.document["thermodynamics"]["gamma"], 1.4)
            self.assertEqual(resolved.source_paths, (case_path.resolve(),))
            self.assertEqual(resolved.extension_paths, {})

    def test_explicit_extensions_are_installed_at_registered_targets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / "config"
            config.mkdir()
            case_path = self._write_modular_case(
                root,
                "  multicomponent: config/multicomponent.yaml\n"
                "  thermodynamics: config/thermodynamics.yaml",
            )
            (config / "multicomponent.yaml").write_text(
                "schema_version: 1\nextension: multicomponent\n"
                "config:\n  mode: inviscid_euler\n  species: [a, b]\n",
                encoding="utf-8",
            )
            (config / "thermodynamics.yaml").write_text(
                "schema_version: 1\nextension: thermodynamics\n"
                "config:\n  model: calorically_perfect\n  gamma: 1.4\n",
                encoding="utf-8",
            )

            resolved = resolve_case_configuration(case_path)

            self.assertNotIn("extensions", resolved.document)
            self.assertEqual(
                resolved.document["physics"]["multicomponent"]["species"],
                ["a", "b"],
            )
            self.assertEqual(resolved.document["thermodynamics"]["gamma"], 1.4)
            self.assertEqual(len(resolved.source_paths), 3)

    def test_duplicate_inline_and_external_configuration_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_path = self._write_modular_case(
                root, "  thermodynamics: thermodynamics.yaml"
            )
            with case_path.open("a", encoding="utf-8") as handle:
                handle.write("thermodynamics:\n  model: calorically_perfect\n")
            (root / "thermodynamics.yaml").write_text(
                "schema_version: 1\nextension: thermodynamics\n"
                "config:\n  model: thermally_perfect\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                CaseConfigurationError, "duplicate configuration for thermodynamics"
            ):
                resolve_case_configuration(case_path)

    def test_extension_cannot_escape_case_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_dir = root / "case"
            case_dir.mkdir()
            case_path = self._write_modular_case(
                case_dir, "  transport: ../transport.yaml"
            )
            (root / "transport.yaml").write_text(
                "schema_version: 1\nextension: transport\nconfig:\n  model: none\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(CaseConfigurationError, "escapes"):
                resolve_case_configuration(case_path)

    def test_sidecar_identity_must_match_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_path = self._write_modular_case(
                root, "  transport: transport.yaml"
            )
            (root / "transport.yaml").write_text(
                "schema_version: 1\nextension: chemistry\nconfig:\n  model: none\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                CaseConfigurationError, "must declare extension: transport"
            ):
                resolve_case_configuration(case_path)

    def test_resolved_case_is_standalone_and_reloadable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_path = self._write_modular_case(
                root, "  chemistry: chemistry.yaml"
            )
            (root / "chemistry.yaml").write_text(
                "schema_version: 1\nextension: chemistry\nconfig:\n  model: none\n",
                encoding="utf-8",
            )
            configuration = resolve_case_configuration(case_path)

            output = write_resolved_case(configuration)
            reloaded = load_yaml(output)

            self.assertEqual(output, root / "resolved_case.yaml")
            self.assertNotIn("extensions", reloaded)
            self.assertEqual(reloaded["chemistry"]["model"], "none")


if __name__ == "__main__":
    unittest.main()
