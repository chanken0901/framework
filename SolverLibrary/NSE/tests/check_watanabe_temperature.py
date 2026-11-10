"""Single-rank Watanabe regression: temperature and SSPRK3 energy budget.

Run with --exe /path/to/nse_cuda --work-dir /new/test/directory.
Requires NumPy. Does not modify any production case or output.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from slf_to_paraview_merged_cropghost import read_slf


def state_sample(q, time, n):
    if q.shape == (n+6, n+6, n+6, 5):
        q = q[3:-3, 3:-3, 3:-3, :]
    if q.shape != (n, n, n, 5) or not np.all(np.isfinite(q)):
        raise ValueError("Invalid state shape or nonfinite SLF")
    rho = q[..., 0]
    if np.any(rho <= 0):
        raise ValueError("Nonpositive density")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        kinetic = .5*np.sum(q[..., 1:4]**2, axis=-1)/rho
        pressure = .4*(q[..., 4]-kinetic)
        if np.any(pressure <= 0):
            raise ValueError("Nonpositive pressure")
        sample = dict(time=float(time), volume_temperature=float(np.mean(pressure/rho)),
                      mass_temperature=float(np.mean(pressure)/np.mean(rho)),
                      mean_density=float(np.mean(rho)), energy=float(np.mean(q[..., 4])),
                      kinetic_energy=float(np.mean(kinetic)))
    if not all(np.isfinite(v) for v in sample.values()):
        raise ValueError("Nonfinite derived statistic or time")
    return sample


def check_samples(samples, end_time):
    samples.sort(key=lambda row: row["time"])
    times = np.array([s["time"] for s in samples])
    if len(times) != 11 or not np.all(np.isfinite(times)) or not np.allclose(
            times, np.linspace(0, end_time, 11), rtol=0, atol=1e-10):
        raise ValueError("Missing, duplicate or unexpected output times (expected 11 snapshots)")
    return {key: max(abs(row[key]/samples[0][key]-1) for row in samples)
            for key in ("volume_temperature", "mass_temperature", "mean_density")}


def check_budget(log, steps, dt, samples):
    coefficients, thermal = [], []
    for line in log.splitlines():
        parts = line.split()
        if parts[:5] == ["#", "forcing", "variance", "dissipation", "thermal_source"]:
            thermal.append([float(v.replace("D", "E")) for v in parts[5:]])
        elif len(parts) > 2 and parts[:2] == ["#", "forcing"] and parts[2].isdigit():
            coefficients.append([float(v.replace("D", "E")) for v in parts[2:]])
    c, h = np.asarray(coefficients), np.asarray(thermal)
    if c.shape != (3*steps, 8) or h.shape != (3*steps, 3):
        raise ValueError("Missing SSPRK3 forcing diagnostics or unexpected retries/substeps")
    if not np.all(np.isfinite(c)) or not np.all(np.isfinite(h)):
        raise ValueError("Nonfinite forcing diagnostic")
    if not np.array_equal(c[:, 0], np.arange(1, 3*steps+1)):
        raise ValueError("Nonsequential forcing evaluations")
    # Columns: evaluation, A, c_d, S, denominator_d, D, epsilon_s, epsilon_d.
    # CPU logs use ES13.5; tolerance includes their rounding, as well as CUDA logs.
    if np.any(h[:, 0] <= 0) or np.any(h[:, 1] < 0) or np.any(c[:, 3] <= 0):
        raise ValueError("Invalid variance, dissipation or solenoidal denominator")
    if np.any(c[:, 2] != 0) or np.any(c[:, 7] != 0):
        raise ValueError("Unexpected dilatational forcing")
    if not np.allclose(h[:, 2], c[:, 5]-h[:, 1], rtol=3e-5, atol=1e-10):
        raise ValueError("Thermal source does not balance pressure dilatation and dissipation")
    expected_a = (h[:, 0]/.04)**(-4)*(.012-c[:, 5])/c[:, 3]
    if not np.allclose(c[:, 1], expected_a, rtol=5e-5, atol=1e-10):
        raise ValueError("Watanabe feedback coefficient mismatch")
    work = c[:, 1]*c[:, 3]
    # SSPRK3 stages are evaluated at q0, q1, q2, with these final weights.
    weights = np.array([1/6, 1/6, 2/3])
    increments = dt*((work+h[:, 2]).reshape(steps, 3) @ weights)
    cumulative = np.concatenate(([0.0], np.cumsum(increments)))
    indices = np.rint(np.array([s["time"] for s in samples])/dt).astype(int)
    observed = np.array([s["energy"]-samples[0]["energy"] for s in samples])
    residual = float(np.max(np.abs(observed-cumulative[indices])))
    scale = dt*float(np.sum((np.abs(work)+np.abs(h[:, 2])).reshape(steps, 3) @ weights))
    if residual > 1e-10+5e-5*scale:
        raise ValueError(f"Total-energy work/heat budget mismatch: {residual:g}")
    return dict(max_energy_budget_residual=residual,
                integrated_work=float(dt*np.sum(work.reshape(steps, 3) @ weights)),
                integrated_heat=float(dt*np.sum(h[:, 2].reshape(steps, 3) @ weights)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    args.exe = args.exe.resolve(strict=True)
    args.work_dir = args.work_dir.resolve()
    args.work_dir.mkdir(parents=True, exist_ok=False)
    reports = []
    for n, dt in ((16, 0.001), (16, 0.0005), (32, 0.001)):
        run = args.work_dir / f"n{n}_dt{dt}"
        run.mkdir()
        steps = round(0.1 / dt)
        source = (Path(__file__).parent / "input_tgv_cuda_forcing_small.dat").read_text()
        changes = {
            'nx = 8': f'nx = {n}', 'ny = 8': f'ny = {n}', 'nz = 8': f'nz = {n}',
            'dt = 1.0e-5': f'dt = {dt}', 't_max = 1.0e-5': 't_max = 0.1',
            'nsteps = 1': f'nsteps = {steps}', 'output_frequency = 0': f'output_frequency = {steps//10}',
            'write_initial = .false.': 'write_initial = .true.',
            'write_meta = .false.': 'write_meta = .true.',
            'mach = 0.1': 'mach = 0.4', 'reynolds = 100.0': 'reynolds = 10.0',
            '"petersen_livescu"': '"watanabe2021"', '"low_wavenumber"': '"full_spectrum"',
            'forcing_fft_backend = "cufft"': 'forcing_fft_backend = "auto"',
            'forcing_target_dissipation = 0.01': 'forcing_target_dissipation = 0.012\n'
                '  forcing_target_velocity_variance = 0.04\n  forcing_feedback_exponent = 4.0',
            'forcing_report_interval = 1': 'forcing_report_interval = 1',
        }
        for old, new in changes.items():
            if old not in source:
                raise RuntimeError(f"Template changed: {old}")
            source = source.replace(old, new)
        (run / "input.dat").write_text(source)
        with (run / "run.log").open("w", encoding="utf-8") as log_file:
            result = subprocess.run([str(args.exe), str(run/'input.dat')],
                                    cwd=run, stdout=log_file, stderr=subprocess.STDOUT, timeout=300)
        if result.returncode:
            raise RuntimeError(f"Solver failed; see {run/'run.log'}")
        samples = []
        for path in sorted(run.rglob("*.slf")):
            slf = read_slf(path)
            samples.append(state_sample(slf.data, slf.time, n))
        drift = check_samples(samples, .1)
        budget = check_budget((run / "run.log").read_text(encoding="utf-8"), steps, dt, samples)
        reports.append(dict(n=n, dt=dt, drift=drift, budget=budget, samples=samples))
        (args.work_dir / "temperature_report.json").write_text(
            json.dumps(reports, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps([{k:r[k] for k in ("n","dt","drift")} for r in reports],indent=2))
    # This is a finite-resolution tolerance, not an assertion of pointwise isothermality.
    if any(r["drift"]["mass_temperature"]>1e-6 or r["drift"]["mean_density"]>1e-12
           or r["drift"]["volume_temperature"]>1e-5 for r in reports):
        raise RuntimeError("Mean internal-temperature or mass drift exceeds tolerance")
    coarse = reports[0]["drift"]["mass_temperature"]
    if any(r["drift"]["mass_temperature"] > coarse+1e-10 for r in reports[1:]):
        raise RuntimeError("Time/grid refinement worsens mass-temperature drift")
    print("Watanabe temperature regression passed")


if __name__ == "__main__":
    main()
