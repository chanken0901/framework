"""Mapped SLF restart matrix; creates isolated test runs, never modifies research runs."""
import argparse
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from slf_to_paraview_merged_cropghost import read_slf


def assemble(directory, step):
    result = np.empty((12, 12, 12, 5))
    covered = np.zeros((12, 12), dtype=bool)
    time = None
    for path in sorted(directory.glob(f"field_{step:06d}*.slf")):
        slf = read_slf(path, allow_nonuniform=True)
        with path.open("rb") as f:
            f.seek(288 + 8*int(np.prod(slf.shape)) + 80)
            mark, *values = struct.unpack("<8s15i15d8s", f.read(196))
        assert mark == b"NSEGRID1" and values[-1] == b"NSREND1\0"
        h = values[:15]
        assert h[1] == step
        ys = slice(h[5]-1, h[6]); zs = slice(h[7]-1, h[8]); g = int(slf.meta[5])
        assert not covered[ys, zs].any()
        covered[ys, zs] = True
        result[:, ys, zs, :] = slf.data[g:g+12, g:g+h[6]-h[5]+1, g:g+h[8]-h[7]+1, :]
        assert time is None or time == slf.time
        time = slf.time
    assert covered.all()
    return result, time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu", type=Path, required=True)
    parser.add_argument("--cuda", type=Path)
    parser.add_argument("--mpi", type=Path)
    parser.add_argument("--mpi-cuda", type=Path)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    commands = {"cpu": [str(args.cpu.resolve())]}
    if args.cuda:
        commands["cuda"] = [str(args.cuda.resolve())]
    for name, exe in (("mpi", args.mpi), ("mpi_cuda", args.mpi_cuda)):
        if exe:
            assert args.mpiexec
            commands[name] = [str(args.mpiexec.resolve()), "-n", "4", str(exe.resolve())]
    work = Path(tempfile.mkdtemp(prefix="nonuniform-restart-", dir=args.work_dir.resolve()))
    print("Artifacts:", work, flush=True)
    env = dict(os.environ, OMP_NUM_THREADS="2", NSE_CUDA_DEVICE_POLICY="fixed")
    template = Path(__file__).with_name("input_nonuniform_production.dat").read_text()
    template = template.replace("output_format='vtr'", "output_format='slf'")
    template = template.replace("nsteps=20", "nsteps=4").replace("output_frequency=20", "output_frequency=2")
    template = template.replace("convective_scheme='keep2'", "convective_scheme='hybrid'")
    template = template.replace("hybrid_smooth_scheme='keep2'", "hybrid_smooth_scheme='keep6'")
    template = template.replace("viscous_scheme='fv2'", "viscous_scheme='central6'")

    def run(name, command, text, error=None):
        folder = work / name
        folder.mkdir()
        inp = folder / "input.dat"
        inp.write_text(text)
        p = subprocess.run(command + [str(inp)], cwd=folder, env=env, capture_output=True, text=True, timeout=90)
        message = p.stdout+p.stderr
        (folder / "run.log").write_text(message)
        if error is None:
            assert p.returncode == 0, message
        else:
            assert p.returncode != 0 and error in message, message
            assert not (folder / "nonuniform-output").exists(), "failed restart wrote output"
        return folder, message

    for boundary in ("periodic", "non_reflecting"):
        text = template
        if boundary != "periodic":
            text = text.replace("boundary_condition='periodic'", """
 boundary_x_min='non_reflecting', boundary_x_max='non_reflecting'
 boundary_y_min='periodic', boundary_y_max='periodic'
 boundary_z_min='periodic', boundary_z_max='periodic'
 boundary_x_min_reference_rho=1, boundary_x_max_reference_rho=1
 boundary_x_min_reference_p=0.7142857142857143, boundary_x_max_reference_p=0.7142857142857143""")
        sources = {}
        for name, cmd in commands.items():
            folder, _ = run(f"{boundary}-source-{name}", cmd, text)
            sources[name] = folder / "nonuniform-output"
        for source, directory in sources.items():
            original, time = assemble(directory, 4)
            snapshot = next(directory.glob("field_000002*.slf"))
            restart = text.replace("&simulation", f"&simulation\n restart_file='{snapshot.as_posix()}'")
            restart = restart.replace("reynolds=100.0", "reynolds=999.0")
            for destination, cmd in commands.items():
                folder, log = run(f"{boundary}-{source}-to-{destination}", cmd, restart)
                resumed, resumed_time = assemble(folder / "nonuniform-output", 4)
                np.testing.assert_allclose(resumed, original, rtol=0, atol=2e-11)
                assert resumed_time == time and "Restart loaded: step=2" in log
                assert not list((folder / "nonuniform-output").glob("field_000000*"))
                print("PASS", boundary, source, "->", destination, flush=True)
        if boundary == "periodic":
            snapshot = next(sources["cpu"].glob("field_000002*.slf"))
            restart = text.replace("&simulation", f"&simulation\n restart_file='{snapshot.as_posix()}'")
            for name, altered, error in (
                ("stretch", restart.replace("0.4,1.0,1.5", "0.5,1.0,1.5"), "bounds/stretch mismatch"),
                ("norm", restart.replace("hybrid_smooth_scheme='keep6'", "hybrid_smooth_scheme='keep2'")
                 .replace("viscous_scheme='central6'", "viscous_scheme='fv2'"), "integration norm mismatch"),
                ("end-step", restart.replace("nsteps=4", "nsteps=2"), "greater than saved"),
                ("end-time", restart.replace("t_max=10.0", "t_max=0.002"), "greater than saved"),
            ):
                run("reject-"+name, commands["cpu"], altered, error)
            broken = work / "broken.slf"
            shutil.copyfile(snapshot, broken)
            with broken.open("r+b") as f:
                f.truncate(broken.stat().st_size-1)
            run("reject-truncated", commands["cpu"], restart.replace(snapshot.as_posix(), broken.as_posix()),
                "size mismatch")
            invalid = work / "invalid-state.slf"
            shutil.copyfile(snapshot, invalid)
            slf = read_slf(snapshot, allow_nonuniform=True)
            g = int(slf.meta[5])
            with invalid.open("r+b") as f:
                f.seek(288+8*(g+slf.shape[0]*(g+slf.shape[1]*g)))
                f.write(struct.pack("<d", float("nan")))
            run("reject-state", commands["cpu"], restart.replace(snapshot.as_posix(), invalid.as_posix()),
                "Nonfinite checkpoint state")
            if "mpi" in sources:
                incomplete = work / "incomplete"
                incomplete.mkdir()
                for src in sources["mpi"].glob("field_000002*.slf"):
                    if "rank00003" not in src.name:
                        shutil.copyfile(src, incomplete / src.name)
                missing = incomplete / "field_000002_rank00000.slf"
                run("reject-rank", commands["cpu"], restart.replace(snapshot.as_posix(), missing.as_posix()),
                    "Missing NSE checkpoint piece")
            repeated = work / "periodic-cpu-to-cpu"
            before = (repeated / "nonuniform-output/field_000004_rank00000.slf").read_bytes()
            p = subprocess.run(commands["cpu"]+[str(repeated / "input.dat")], cwd=repeated,
                               env=env, capture_output=True, text=True, timeout=30)
            assert p.returncode != 0 and "new output directory" in p.stdout+p.stderr
            assert (repeated / "nonuniform-output/field_000004_rank00000.slf").read_bytes() == before
            try:
                read_slf(snapshot)
            except ValueError as exc:
                assert "Nonuniform SLF" in str(exc)
            else:
                raise AssertionError("uniform-grid postprocessing must reject mapped data")
    print("All mapped SLF restart checks passed.", flush=True)


if __name__ == "__main__":
    main()
