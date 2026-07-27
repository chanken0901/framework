from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from global_case_index import (  # noqa: E402
    rebuild_case_index,
    sync_case_document,
)


def _case(case_id: str, model: str, profile: str, nx: int) -> dict:
    return {
        "case_id": case_id,
        "case_label": f"{model}_{case_id}",
        "status": "planned",
        "description": "test",
        "physics": {
            "model": model,
            model: {
                "alpha": 0.05,
            },
        },
        "flow": {"type": "quantum_taylor_green"},
        "grid": {"nx": nx, "ny": nx, "nz": nx},
        "time": {"dt": 1.0e-4, "nsteps": 100},
        "solver": {"profile": profile, "processes": 1},
    }


class GlobalCaseIndexTests(unittest.TestCase):
    def test_sync_lists_multiple_cases_and_updates_conditions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            index = root / "case_index.csv"
            case1_path = root / "gpe_case0001" / "cases" / "case0001" / "case.yaml"
            case2_path = root / "gpe_case0002" / "cases" / "case0002" / "case.yaml"

            sync_case_document(
                _case("case0001", "gpe", "cuda_single", 64),
                index_path=index,
                model="gpe",
                profile="cuda_single",
                processes=1,
                environment="gpe_case0001",
                case_path=case1_path,
            )
            sync_case_document(
                _case("case0002", "gpe", "cpu_mpi_fftw", 256),
                index_path=index,
                model="gpe",
                profile="cpu_mpi_fftw",
                processes=8,
                environment="gpe_case0002",
                case_path=case2_path,
            )

            with index.open("r", newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["case_id"] for row in rows], ["case0001", "case0002"])
            self.assertEqual(rows[1]["grid.nx"], "256")
            self.assertIn("physics.gpe.alpha", rows[1])

            registered = rows[1]["registered_at_utc"]
            sync_case_document(
                _case("case0002", "gpe", "cpu_mpi_fftw", 512),
                index_path=index,
                model="gpe",
                profile="cpu_mpi_fftw",
                processes=8,
                environment="gpe_case0002",
                case_path=case2_path,
            )
            with index.open("r", newline="", encoding="utf-8") as handle:
                updated = list(csv.DictReader(handle))
            self.assertEqual(len(updated), 2)
            self.assertEqual(updated[1]["grid.nx"], "512")
            self.assertEqual(updated[1]["registered_at_utc"], registered)

    def test_rebuild_scans_generated_environments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            environment = root / "gpe_case0001"
            case_path = environment / "cases" / "case0001" / "case.yaml"
            case_path.parent.mkdir(parents=True)
            case_path.write_text(
                json.dumps(_case("case0001", "gpe", "cuda_single", 64)),
                encoding="utf-8",
            )
            (environment / "environment.lock.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "model": "gpe",
                        "profile": "cuda_single",
                        "processes": 1,
                        "case_id": "case0001",
                        "case_directory": "cases/case0001",
                    }
                ),
                encoding="utf-8",
            )

            count = rebuild_case_index(root, root / "case_index.csv")

            self.assertEqual(count, 1)
            with (root / "case_index.csv").open(
                "r", newline="", encoding="utf-8"
            ) as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["case_key"], "gpe:case0001")
            self.assertEqual(row["environment"], "gpe_case0001")


if __name__ == "__main__":
    unittest.main()
