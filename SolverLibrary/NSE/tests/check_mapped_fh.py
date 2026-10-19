"""End-to-end mapped LLNS comparison and restart matrix (isolated scratch runs)."""
import argparse
import os
from pathlib import Path
import subprocess
import shutil
import tempfile
import numpy as np
from check_nonuniform_restart import assemble


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('cpu', 'cuda', 'mpi', 'mpi-cuda', 'mpiexec'):
        parser.add_argument('--'+option, type=Path, required=option == 'cpu')
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--boundary-case', choices=('periodic','no_external_noise','thermal_reservoir','mixed'), default='periodic')
    parser.add_argument('--open-axis', choices=('x','y','z','all'), default='all')
    parser.add_argument('--uniform', action='store_true')
    args = parser.parse_args()
    mpi_launcher = str(args.mpiexec.resolve()) if args.mpiexec else shutil.which('mpiexec')
    if (args.mpi or args.mpi_cuda) and not mpi_launcher:
        parser.error('MPI tests require --mpiexec or mpiexec on PATH')
    commands = {'cpu': [str(args.cpu.resolve())]}
    if args.cuda:
        commands['cuda'] = [str(args.cuda.resolve())]
    for name, exe in (('mpi', args.mpi), ('mpi_cuda', args.mpi_cuda)):
        if exe:
            commands[name] = [mpi_launcher, '-n', '4', str(exe.resolve())]
    root = Path(tempfile.mkdtemp(prefix='mapped-fh-', dir=args.work_dir.resolve()))
    print('Artifacts:', root, flush=True)
    text = Path(__file__).with_name('input_nonuniform_production.dat').read_text()
    text = text.replace("output_format='vtr'", "output_format='slf'")
    text = text.replace('nsteps=20', 'nsteps=4').replace('output_frequency=20', 'output_frequency=2')
    text = text.replace("viscous_scheme='fv2'", "viscous_scheme='central6'")
    text = text.replace('&nse', '&nse\n fh_enabled=.true.,fh_boltzmann_number=1.e-8,fh_seed=13579')
    if args.uniform:
        text = text.replace("grid_mapping='sinh', grid_stretch=0.4,1.0,1.5", "grid_mapping='uniform'")
    if args.boundary_case != 'periodic':
        faces = ('x_min','x_max','y_min','y_max','z_min','z_max')
        lines = []
        modes = []
        for index, face in enumerate(faces):
            opened = args.open_axis == 'all' or face.startswith(args.open_axis)
            lines.append("boundary_%s='%s'" % (face, 'non_reflecting' if opened else 'periodic'))
            modes.append(int(opened and (args.boundary_case == 'thermal_reservoir' or
                                        (args.boundary_case == 'mixed' and index % 2 == 1))))
            if opened:
                for key, value in (('rho',1),('u',0),('v',0),('w',0),('p',1)):
                    lines.append('boundary_%s_reference_%s=%s' % (face,key,value))
        lines.append('fh_boundary_model='+','.join(map(str,modes)))
        text = text.replace("boundary_condition='periodic'", '\n'.join(lines))
    env = dict(os.environ, OMP_NUM_THREADS='2', NSE_CUDA_DEVICE_POLICY='fixed')

    def run(name, cmd, content):
        folder = root/name
        folder.mkdir()
        inp = folder/'input.dat'; inp.write_text(content)
        p = subprocess.run(cmd+[str(inp)], cwd=folder, env=env, capture_output=True, text=True, timeout=90)
        (folder/'run.log').write_text(p.stdout+p.stderr)
        assert p.returncode == 0, p.stdout+p.stderr
        return folder/'nonuniform-output', p.stdout

    for scheme in ('keep2', 'keep6', 'weno5z_roe', 'hybrid2', 'hybrid6'):
        for adaptive in (False, True):
            content = text.replace("convective_scheme='keep2'", "convective_scheme='%s'" % ('hybrid' if scheme.startswith('hybrid') else scheme))
            if scheme == 'hybrid6':
                content = content.replace("hybrid_smooth_scheme='keep2'", "hybrid_smooth_scheme='keep6'")
            if adaptive:
                content = content.replace('use_fixed_dt=.true.', 'use_fixed_dt=.false.')
            label = scheme+('-adaptive' if adaptive else '-fixed')
            reference = None
            source = None
            for backend, cmd in commands.items():
                output, _ = run(label+'-'+backend, cmd, content)
                state, time = assemble(output, 4)
                if reference is None:
                    reference, ref_time, source = state, time, output
                else:
                    np.testing.assert_allclose(state, reference, atol=3e-11, rtol=1e-11)
                    np.testing.assert_allclose(time, ref_time, atol=1e-13, rtol=1e-12)
                restart = content.replace('&simulation', "&simulation\n restart_file='%s'" % (source/'field_000002_rank00000.slf').as_posix())
                resumed, log = run(label+'-restart-'+backend, cmd, restart)
                state, time = assemble(resumed, 4)
                np.testing.assert_allclose(state, reference, atol=3e-11, rtol=1e-11)
                np.testing.assert_allclose(time, ref_time, atol=1e-13, rtol=1e-12)
                assert 'Restart loaded: step=2' in log
                print('PASS', label, backend, 'and restart', flush=True)
    print('Mapped LLNS backend/restart matrix passed', flush=True)


if __name__ == '__main__':
    main()
