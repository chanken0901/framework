from __future__ import annotations

import json
import struct
import sys
import tempfile
import os
import subprocess
import unittest
from pathlib import Path

import numpy as np


TOOLS_DIR = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from nse_prepare_imported_turbulence import (  # noqa: E402
    prepare_imported_turbulence,
    read_parameters, parameters_from_input, write_global_slf, PARAMETER_FORMAT,
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
        stream.write(PARAMETER_FORMAT.pack(b"NSEPAR1\0",1.4,1767.7669529663692,.72,1.,.3,
                                           b"central6".ljust(32,b" ")))


class PrepareImportedTurbulenceTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get("NSE_TEST_CPU"), "set NSE_TEST_CPU for solver integration")
    def test_solver_restart_and_import_inherit_physics(self):
        commands = {"cpu": [os.environ["NSE_TEST_CPU"]]}
        if os.environ.get("NSE_TEST_CUDA"):
            commands["cuda"] = [os.environ["NSE_TEST_CUDA"]]
        for name in ("MPI", "MPI_CUDA"):
            if os.environ.get("NSE_TEST_"+name):
                commands[name.lower()] = [os.environ["NSE_TEST_MPIEXEC"], "-n", "4", os.environ["NSE_TEST_"+name]]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/"source"; source.mkdir()
            text = '''&simulation
 equation='NSE', initial_condition='taylor_green', nx=8,ny=8,nz=8,nghost=3,
 x_max=6.283185307179586,y_max=6.283185307179586,z_max=6.283185307179586,
 dt=.0001,t_max=1,nsteps=2,use_fixed_dt=.true.,output_frequency=1,
 output_dir='out',write_initial=.true.,write_meta=.true.
/
&nse
 gamma=1.4,reynolds=1767.7669529663692,prandtl=.72,rho0=1,mach=.3,
 convective_scheme='hybrid',hybrid_smooth_scheme='keep6',hybrid_shock_scheme='weno5z_roe',
 viscous_scheme='central6',boundary_condition='periodic'
/
'''
            env = dict(os.environ, OMP_NUM_THREADS="2", NSE_CUDA_DEVICE_POLICY="fixed")
            def run(directory, input_text, command, success=True):
                (directory/"input.dat").write_text(input_text)
                result = subprocess.run(command+[str(directory/"input.dat")],cwd=directory,env=env,
                                        capture_output=True,text=True,timeout=60)
                self.assertEqual(result.returncode == 0,success,result.stdout+result.stderr)
                return result.stdout+result.stderr
            run(source,text,commands["cpu"])
            portable = root/"source.slf"
            prepare_imported_turbulence(source/"out",portable,step="1")
            reference = root/"reference.slf"
            prepare_imported_turbulence(source/"out",reference,step="2")
            expected = read_slf(reference).data
            for name, command in commands.items():
                for restart in (True,False):
                    target = root/f"{name}-{restart}"; target.mkdir()
                    bad = text.replace("gamma=1.4,reynolds=1767.7669529663692,prandtl=.72,rho0=1,mach=.3",
                                       "gamma=1.6,reynolds=100,prandtl=.9,rho0=2,mach=.7")
                    if restart:
                        bad = bad.replace("&simulation",f"&simulation\n restart_file='{portable.as_posix()}',")
                    else:
                        bad = bad.replace("initial_condition='taylor_green'", "initial_condition='imported_turbulence'")
                        bad = bad.replace("nsteps=2", "nsteps=1")
                        bad = bad.replace("&nse",f"&nse\n imported_turbulence_file='{portable.as_posix()}', "
                                          "imported_turbulence_mode='tile',imported_turbulence_blend_cells=0,")
                    log = run(target,bad,command)
                    self.assertIn("# inherited gamma,Re,Pr,rho0,mach:",log)
                    merged = root/f"{name}-{restart}.slf"
                    prepare_imported_turbulence(target/"out",merged)
                    self.assertEqual(read_parameters(merged),read_parameters(portable))
                    np.testing.assert_allclose(read_slf(merged).data,expected,rtol=0,atol=2e-11)
            # Direct namelist invocation must not silently use destination Re.
            legacy = root/"legacy.slf"
            legacy.write_bytes(portable.read_bytes()[:-80])
            rejected = root/"rejected"; rejected.mkdir()
            log = run(rejected,text.replace("&simulation",f"&simulation\n restart_file='{legacy.as_posix()}',"),
                      commands["cpu"],False)
            self.assertIn("Source SLF lacks NSE parameters",log)

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
            self.assertEqual(read_parameters(output), read_parameters(source / "field_000007_rank00000.slf"))
            rank1 = source/"field_000007_rank00001.slf"
            rank1.write_bytes(rank1.read_bytes()[:-80] + PARAMETER_FORMAT.pack(
                b"NSEPAR1\0",1.4,100.,.72,1.,.3,b"central6".ljust(32,b" ")))
            with self.assertRaisesRegex(ValueError,"rank files disagree"):
                prepare_imported_turbulence(source,root/"bad.slf")
            self.assertFalse((root/"bad.slf").exists())

    def test_legacy_requires_explicit_resolved_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, target = root/"field_000007.slf", root/"portable.slf"
            field = np.zeros((2,2,2,5)); field[:,:,:,0]=1; field[:,:,:,4]=2.5
            write_global_slf(source, field, list(NAMES), 7, .25, (0,0,0), (1,1,1))
            with self.assertRaisesRegex(ValueError, "source SLF lacks"):
                prepare_imported_turbulence(source,target)
            self.assertFalse(target.exists())
            inp = root/"original.dat"
            inp.write_text('&nse\n gamma=1.4, reynolds=1767.7669529663692, prandtl=.72, rho0=1, mach=.3,\n viscous_scheme="central6"\n/\n')
            prepare_imported_turbulence(source,target,source_input=inp)
            self.assertEqual(read_parameters(target), parameters_from_input(inp))
            with self.assertRaisesRegex(ValueError, "gamma disagrees"):
                prepare_imported_turbulence(target,root/"other.slf",gamma=1.6)


if __name__ == "__main__":
    unittest.main()
