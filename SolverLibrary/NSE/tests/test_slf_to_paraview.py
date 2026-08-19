from __future__ import annotations

import json
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


SOLVER_LIBRARY = Path(__file__).resolve().parents[2]
CONVERTERS = (
    SOLVER_LIBRARY / "NSE" / "tools" / "slf_to_paraview_merged_cropghost.py",
    SOLVER_LIBRARY / "GPE" / "gp3d" / "tools" / "slf_to_paraview_merged_cropghost.py",
)


def write_slf(path: Path, step: int) -> None:
    shape = (2, 2, 2, 2)
    meta = np.zeros(8, dtype="<i4")
    meta[0] = step
    data = np.zeros(shape, dtype="<f8", order="F")
    data[..., 0] = 1.0
    data[..., 1] = 0.5
    with path.open("wb") as handle:
        handle.write(b"SLF1\0\0\0\0")
        handle.write(struct.pack("<iii", 1, 2, 3))
        handle.write(np.asarray(shape, dtype="<i4").tobytes())
        handle.write(meta.tobytes())
        handle.write(struct.pack("<d", float(step)))
        handle.write(struct.pack("<6d", 0.0, 2.0, 0.0, 2.0, 0.0, 2.0))
        handle.write(struct.pack("<i", 2))
        for name in ("psi_real", "psi_imag"):
            handle.write(name.encode("ascii").ljust(32, b"\0"))
        handle.write(data.tobytes(order="F"))


class SlfToParaViewTests(unittest.TestCase):
    def test_both_converters_support_steps_and_inspect_only(self) -> None:
        for converter in CONVERTERS:
            with self.subTest(converter=converter):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    output = root / "output"
                    output.mkdir()
                    write_slf(output / "field_000000.slf", 0)
                    write_slf(output / "field_000100.slf", 100)
                    (output / "meta.json").write_text(
                        json.dumps(
                            {
                                "equation": "GPE",
                                "grid": [2, 2, 2],
                                "parallel": {"mpi_nprocs": 1},
                            }
                        ),
                        encoding="utf-8",
                    )
                    destination = root / "paraview"
                    completed = subprocess.run(
                        [
                            sys.executable,
                            str(converter),
                            str(output),
                            "--output-dir",
                            str(destination),
                            "--derive",
                            "gpe",
                            "--steps",
                            "latest",
                            "--fields",
                            "density,phase",
                            "--inspect-only",
                        ],
                        check=False,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(completed.returncode, 0, completed.stderr)
                    self.assertIn("selected steps: [100]", completed.stdout)
                    self.assertFalse(destination.exists())

    def test_gpe_converter_writes_selected_vti_and_pvd(self) -> None:
        converter = CONVERTERS[1]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            write_slf(output / "field_000000.slf", 0)
            write_slf(output / "field_000100.slf", 100)
            (output / "meta.json").write_text(
                json.dumps(
                    {
                        "equation": "GPE",
                        "grid": [2, 2, 2],
                        "parallel": {"mpi_nprocs": 1},
                    }
                ),
                encoding="utf-8",
            )
            destination = root / "paraview"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(converter),
                    str(output),
                    "--output-dir",
                    str(destination),
                    "--derive",
                    "gpe",
                    "--steps",
                    "100",
                    "--fields",
                    "density,phase",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((destination / "field_000100.vti").is_file())
            self.assertTrue((destination / "collection.pvd").is_file())
            self.assertFalse((destination / "field_000000.vti").exists())

    def test_global_layout_is_complete_even_with_mpi_metadata(self) -> None:
        converter = CONVERTERS[0]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            write_slf(output / "field_000100.slf", 100)
            (output / "meta.json").write_text(
                json.dumps(
                    {
                        "equation": "GPE",
                        "grid": [2, 2, 2],
                        "parallel": {"mpi_enabled": True, "mpi_nprocs": 4},
                    }
                ),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(converter),
                    str(output),
                    "--derive",
                    "gpe",
                    "--layout",
                    "global",
                    "--steps",
                    "latest",
                    "--inspect-only",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("selected steps: [100]", completed.stdout)


if __name__ == "__main__":
    unittest.main()
