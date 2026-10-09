"""Short end-to-end checks of executables, physical VTK output and backend agreement.

Run with --cpu EXE [--cuda EXE] [--mpi EXE --mpiexec EXE].
All generated inputs/outputs are kept in a new directory below --work-dir.
No research case is modified. Only the Python standard library is required.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET


def read_output(directory, step, uniform_flow=False, mapped_keep6=False, boundary="periodic"):
    master = ET.parse(directory / f"field_{step:06d}.pvtr").getroot()
    field = [None] * 12**3
    volumes = [None] * len(field)
    time = None
    for item in master.findall(".//Piece"):
        path = directory / item.attrib["Source"]
        root = ET.parse(path).getroot()
        piece = root.find(".//Piece")
        ext = list(map(int, piece.attrib["Extent"].split()))
        assert ext == list(map(int, item.attrib["Extent"].split()))
        coords = [list(map(float, a.text.split())) for a in piece.findall("Coordinates/DataArray")]
        data = {a.attrib["Name"]: list(map(float, a.text.split()))
                for a in piece.findall("CellData/DataArray")}
        assert ("integration_weight" in data) == mapped_keep6
        values = list(zip(*(data[name] for name in ("rho", "rho_u", "rho_v", "rho_w", "rho_E"))))
        n = 0
        for axis, beta in enumerate((0.4, 1.0, 1.5)):
            assert len(coords[axis]) == ext[2*axis+1] - ext[2*axis] + 1
            for offset, x in enumerate(coords[axis]):
                i = ext[2*axis] + offset
                expected = math.pi * (1 + math.sinh(beta*(2*i/12-1))/math.sinh(beta))
                assert abs(x-expected) < 2e-14, (axis, i, x, expected)
        for k in range(ext[4], ext[5]):
            for j in range(ext[2], ext[3]):
                for i in range(ext[0], ext[1]):
                    idx = i + 12*(j+12*k)
                    assert field[idx] is None, "overlapping output cells"
                    q = values[n]
                    assert all(map(math.isfinite, q)) and q[0] > 0
                    pressure = 0.4*(q[4] - sum(v*v for v in q[1:4])/(2*q[0]))
                    assert pressure > 0, (idx, pressure)
                    field[idx] = q
                    volumes[idx] = math.prod(coords[a][p-ext[2*a]+1]-coords[a][p-ext[2*a]]
                                              for a, p in enumerate((i,j,k)))
                    if mapped_keep6:
                        def metric(axis, p):
                            beta = (0.4, 1.0, 1.5)[axis]
                            edges = [math.pi*(1+math.sinh(beta*(2*m/12-1))/math.sinh(beta))
                                     for m in range(13)]
                            centers = [(a+b)/2 for a,b in zip(edges, edges[1:])]
                            def center(m):
                                if axis != 0 or boundary == "periodic":
                                    return centers[m % 12] + (m//12)*2*math.pi
                                if m < 0:
                                    return -centers[-m-1]
                                if m >= 12:
                                    return 4*math.pi-centers[23-m]
                                return centers[m]
                            return sum(c*(center(p+s)-center(p-s))
                                       for s,c in enumerate((.75,-.15,1/60),1))
                        weight = math.prod(metric(a,p) for a,p in enumerate((i,j,k)))
                        assert weight > 0 and abs(data["integration_weight"][n]-weight) < 2e-14
                        volumes[idx] = weight
                    if uniform_flow:
                        expected_q = (1.2, 1.2*.13, -1.2*.04, 1.2*.02,
                                      .9/.4 + .5*1.2*(.13**2+.04**2+.02**2))
                        assert max(abs(a-b) for a,b in zip(q, expected_q)) < 2e-11, (idx,q)
                    elif step == 0:
                        x, y, z = [0.5*(coords[a][p-ext[2*a]]+coords[a][p-ext[2*a]+1])
                                   for a, p in enumerate((i,j,k))]
                        assert abs(q[1] - .1*math.sin(x)*math.cos(y)*math.cos(z)) < 2e-14
                    n += 1
        assert n == len(values)
        piece_time = float(root.find(".//FieldData/DataArray[@Name='TimeValue']").text)
        assert int(root.find(".//FieldData/DataArray[@Name='Step']").text) == step
        assert time is None or time == piece_time
        time = piece_time
    assert None not in field and None not in volumes
    return field, volumes, time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu", required=True, type=Path)
    parser.add_argument("--cuda", type=Path)
    parser.add_argument("--mpi", type=Path)
    parser.add_argument("--mpi-cuda", type=Path)
    parser.add_argument("--shared-gpu", action="store_true",
                        help="Test four MPI ranks on one GPU; not a multi-GPU performance test")
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--uniform-flow", action="store_true",
                        help="Check moving free-stream preservation, including inlet/outlet boundaries")
    args = parser.parse_args()
    if args.steps < 1:
        parser.error("--steps must be positive")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="nonuniform-production-", dir=args.work_dir.resolve()))
    print("Artifacts:", work, flush=True)
    template = Path(__file__).with_name("input_nonuniform_production.dat").read_text()
    if args.uniform_flow:
        template = template.replace("initial_condition='taylor_green'", "initial_condition='uniform_flow'")
        template = template.replace("&nse", "&nse\n uniform_state=1.2,0.13,-0.04,0.02,0.9")
    template = template.replace("nsteps=20", f"nsteps={args.steps}").replace("t_max=10.0", "t_max=1000.0")
    # The final state must be written even away from a regular output interval.
    template = template.replace("output_frequency=20", f"output_frequency={args.steps+1}")
    commands = {"cpu": [str(args.cpu.resolve())]}
    if args.cuda:
        commands["cuda"] = [str(args.cuda.resolve())]
    if args.mpi:
        if not args.mpiexec:
            parser.error("--mpi requires --mpiexec")
        commands["mpi4"] = [str(args.mpiexec.resolve()), "-n", "4", str(args.mpi.resolve())]
    if args.mpi_cuda:
        if not args.mpiexec:
            parser.error("--mpi-cuda requires --mpiexec")
        commands["mpi_cuda4"] = [str(args.mpiexec.resolve()), "-n", "4", str(args.mpi_cuda.resolve())]
    env = dict(os.environ, OMP_NUM_THREADS="2")
    if args.shared_gpu:
        env["NSE_CUDA_DEVICE_POLICY"] = "fixed"
    for scheme in ("keep2", "keep6", "weno5z_roe", "hybrid", "hybrid6"):
        mapped_norm = scheme in {"keep6", "hybrid6"}
        for viscous in (("none", "fv2", "central6") if mapped_norm else ("none", "fv2")):
            for boundary in ("periodic", "mixed"):
                reference = None
                for backend, command in commands.items():
                    run = work / f"{scheme}-{viscous}-{boundary}-{backend}"
                    run.mkdir()
                    text = template.replace("convective_scheme='keep2'", f"convective_scheme='{scheme}'")
                    if scheme == "hybrid6":
                        text = text.replace("convective_scheme='hybrid6'", "convective_scheme='hybrid'")
                        text = text.replace("hybrid_smooth_scheme='keep2'", "hybrid_smooth_scheme='keep6'")
                    text = text.replace("viscous_scheme='fv2'", f"viscous_scheme='{viscous}'")
                    # Exercise adaptive dt as well as fixed dt with CPU/GPU/MPI parity.
                    if boundary == "mixed":
                        text = text.replace("use_fixed_dt=.true.", "use_fixed_dt=.false.")
                        text = text.replace("boundary_condition='periodic'", """
 boundary_x_min='non_reflecting', boundary_x_max='non_reflecting'
 boundary_y_min='periodic', boundary_y_max='periodic'
 boundary_z_min='periodic', boundary_z_max='periodic'
 boundary_x_min_reference_rho=1, boundary_x_max_reference_rho=1
 boundary_x_min_reference_p=0.7142857142857143, boundary_x_max_reference_p=0.7142857142857143""")
                    inp = run / "input.dat"
                    if args.uniform_flow and boundary == "mixed":
                        text = text.replace("boundary_x_min='non_reflecting'", "boundary_x_min='dirichlet'")
                        text = text.replace("reference_rho=1", "reference_rho=1.2")
                        text = text.replace("reference_p=0.7142857142857143", "reference_p=0.9")
                        text = text.replace("boundary_x_min_reference_rho", """
 boundary_x_min_reference_u=0.13, boundary_x_max_reference_u=0.13
 boundary_x_min_reference_v=-0.04, boundary_x_max_reference_v=-0.04
 boundary_x_min_reference_w=0.02, boundary_x_max_reference_w=0.02
 boundary_x_min_reference_rho""")
                    inp.write_text(text)
                    result = subprocess.run(command+[str(inp)], cwd=run, env=env,
                                            text=True, capture_output=True, timeout=180)
                    (run / "run.log").write_text(result.stdout+result.stderr)
                    assert result.returncode == 0, (run, result.stdout[-2000:], result.stderr[-2000:])
                    assert "completed successfully" in result.stdout
                    out = run / "nonuniform-output"
                    meta = json.loads((out / "meta.json").read_text())
                    assert meta["spacing"] is None and meta["format"] == "vtr"
                    if mapped_norm:
                        assert meta["state_representation"] == "mapped_grid_point_values"
                    assert not list(out.glob("*.slf"))
                    initial, volume, _ = read_output(out, 0, args.uniform_flow, mapped_norm, boundary)
                    final, _, time = read_output(out, args.steps, args.uniform_flow, mapped_norm, boundary)
                    assert time > 0
                    if boundary == "periodic":
                        for v in range(5):
                            before = math.fsum(q[v]*dv for q, dv in zip(initial, volume))
                            after = math.fsum(q[v]*dv for q, dv in zip(final, volume))
                            assert abs(after-before) < 2e-11*max(1, abs(before)), (v, before, after)
                    if reference is None:
                        reference = final, time
                    else:
                        error = max(abs(a-b) for q, qr in zip(final, reference[0]) for a,b in zip(q,qr))
                        assert error < 2e-10 and abs(time-reference[1]) < 2e-12, (run, error, time, reference[1])
                    print("PASS", run.name, f"time={time:.8g}", flush=True)
    # A physical end time reached between output intervals also saves the last state.
    for backend, command in commands.items():
        run = work / f"end-time-{backend}"
        run.mkdir()
        inp = run / "input.dat"
        inp.write_text(template.replace("t_max=1000.0", "t_max=0.0035")
                       .replace(f"nsteps={args.steps}", "nsteps=100"))
        result = subprocess.run(command+[str(inp)], cwd=run, env=env,
                                text=True, capture_output=True, timeout=60)
        (run / "run.log").write_text(result.stdout+result.stderr)
        assert result.returncode == 0, result.stderr
        _, _, time = read_output(run / "nonuniform-output", 4, args.uniform_flow)
        assert abs(time-.0035) < 1e-15
        print("PASS end-time", backend, flush=True)
    # Direct namelist callers also receive explicit rejection before any output.
    rejections = [("output_format='vtr'", "output_format='slf'"),
                     ("hybrid_smooth_scheme='keep2'", "hybrid_smooth_scheme='unsupported', convective_scheme='hybrid'"),
                     ("viscous_scheme='fv2'", "viscous_scheme='central6'")]
    if args.uniform_flow:
        rejections += [("uniform_state=1.2,0.13,-0.04,0.02,0.9", state) for state in
                       ("", "uniform_state=-1,0,0,0,1", "uniform_state=1,0,0,0,0")]
    for old, new in rejections:
        run = Path(tempfile.mkdtemp(prefix="reject-", dir=work))
        inp = run / "input.dat"
        inp.write_text(template.replace(old,new))
        result = subprocess.run(commands["cpu"]+[str(inp)], cwd=run, env=env,
                                text=True, capture_output=True, timeout=30)
        message = result.stderr + result.stdout
        assert result.returncode != 0 and any(s in message for s in ("Nonuniform", "Mapped CENTRAL6", "uniform_state"))
        assert not (run / "nonuniform-output").exists()
    print("All nonuniform production checks passed.", flush=True)


if __name__ == "__main__":
    main()
