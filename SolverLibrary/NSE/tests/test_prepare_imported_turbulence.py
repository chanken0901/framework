from __future__ import annotations

import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


TOOLS_DIR = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from nse_prepare_imported_turbulence import (  # noqa: E402
    prepare_imported_turbulence,
)
from slf_to_paraview_merged_cropghost import read_slf  # noqa: E402


NAMES = ("rho", "rho_u", "rho_v", "rho_w", "rho_E")


def _state(global_i: int, global_j: int, global_k: int) -> np.ndarray:
    rho = 1.0
    u = float(global_i + 10 * global_j + 100 * global_k)
    v = 0.0
    w = 0.0
    pressure = 1.0
    return np.asarray(
        [rho, rho * u, rho * v, rho * w, pressure / 0.4 + 0.5 * rho * u * u]
    )


def _write_rank_slf(path: Path, rank: int, j_start: int, j_end: int) -> None:
    nghost = 1
    nx, ny, nz = 4, j_end - j_start + 1, 2
    field = np.full(
        (nx + 2 * nghost, ny + 2 * nghost, nz + 2 * nghost, 5),
        -999.0,
        dtype="<f8",
        order="F",
    )
    for local_k in range(nz):
        for local_j in range(ny):
            for local_i in range(nx):
                field[
                    local_i + nghost,
                    local_j + nghost,
                    local_k + nghost,
                    :,
                ] = _state(local_i + 1, j_start + local_j, local_k + 1)

    shape = np.asarray(field.shape, dtype="<i4")
    metadata = np.asarray([7, rank, 4, 4, 2, nghost, 2, 0], dtype="<i4")
    with path.open("wb") as stream:
        stream.write(b"SLF1\x00\x00\x00\x00")
        stream.write(struct.pack("<iii", 1, 2, 4))
        stream.write(shape.tobytes())
        stream.write(metadata.tobytes())
        stream.write(struct.pack("<d", 0.25))
        stream.write(struct.pack("<6d", 0.0, 4.0, 0.0, 4.0, 0.0, 2.0))
        stream.write(struct.pack("<i", 5))
        for name in NAMES:
            stream.write(name.encode("ascii").ljust(32, b" "))
        stream.write(field.tobytes(order="F"))


class PrepareImportedTurbulenceTests(unittest.TestCase):
    def test_rank_files_are_merged_and_ghost_cells_removed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            _write_rank_slf(source / "field_000007_rank00000.slf", 0, 1, 2)
            _write_rank_slf(source / "field_000007_rank00001.slf", 1, 3, 4)
            (source / "meta.json").write_text(
                json.dumps(
                    {
                        "equation": "NSE",
                        "grid": [4, 4, 2],
                        "origin": [0.0, 0.0, 0.0],
                        "spacing": [1.0, 1.0, 1.0],
                        "parallel": {
                            "mpi_enabled": True,
                            "mpi_nprocs": 2,
                            "decomposition": "x-global_yz-block",
                            "rank_ranges": [
                                {
                                    "rank": 0,
                                    "i_start": 1,
                                    "i_end": 4,
                                    "j_start": 1,
                                    "j_end": 2,
                                    "k_start": 1,
                                    "k_end": 2,
                                },
                                {
                                    "rank": 1,
                                    "i_start": 1,
                                    "i_end": 4,
                                    "j_start": 3,
                                    "j_end": 4,
                                    "k_start": 1,
                                    "k_end": 2,
                                },
                            ],
                        },
                    }
                ),
                encoding="utf-8",
            )
            output = root / "initial" / "turbulence.slf"

            step, shape = prepare_imported_turbulence(source, output)

            self.assertEqual(step, 7)
            self.assertEqual(shape, (4, 4, 2))
            result = read_slf(output)
            self.assertEqual(result.shape, (4, 4, 2, 5))
            self.assertEqual(int(result.meta[5]), 0)
            self.assertEqual(int(result.meta[6]), 1)
            np.testing.assert_allclose(result.data[2, 3, 1, :], _state(3, 4, 2))


if __name__ == "__main__":
    unittest.main()
