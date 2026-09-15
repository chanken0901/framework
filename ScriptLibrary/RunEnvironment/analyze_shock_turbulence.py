#!/usr/bin/env python3
"""Dedicated generated-case entry point for shock/turbulence diagnostics."""
import argparse
import subprocess
import sys
from postprocess_case import _case_context, _resolve, _display_command, PostprocessError


def build_command(args):
    context = _case_context(args)
    if context['model'] != 'nse':
        raise PostprocessError('Shock analysis supports only single-component NSE')
    tool = context['root']/'SolverLibrary/NSE/tools/nse_shock_analysis.py'
    if not tool.is_file():
        raise PostprocessError('Shock analysis tool missing; regenerate/update the run environment')
    output = _resolve(context['root'],args.output_dir,context['case_root']/'shock_analysis')
    command = [sys.executable,str(tool),str(context['input_dir']), '--meta',str(context['meta']),
               '--output-dir',str(output),'--gamma',str(context['gamma']),
               '--direction',str(args.direction),'--steps',args.steps,'--layout',args.layout,
               '--min-pressure-ratio',str(args.min_pressure_ratio)]
    for key in ('search','upstream','downstream'):
        command += ['--'+key,*map(str,getattr(args,key))]
    if args.max_shift is not None:
        command += ['--max-shift',str(args.max_shift)]
    return command


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case-directory'); p.add_argument('--input-dir'); p.add_argument('--meta')
    p.add_argument('--gamma',type=float); p.add_argument('--output-dir')
    p.add_argument('--direction',type=int,choices=(-1,1),required=True)
    for name in ('search','upstream','downstream'):
        p.add_argument('--'+name,nargs=2,type=float,required=True)
    p.add_argument('--steps',default='all')
    p.add_argument('--layout',choices=('auto','global','rank'),default='auto')
    p.add_argument('--max-shift',type=float)
    p.add_argument('--min-pressure-ratio',type=float,default=1.05)
    p.add_argument('--dry-run',action='store_true')
    p.set_defaults(task='shock',reynolds=None)
    args = p.parse_args(argv)
    try:
        command = build_command(args)
        print('[CMD]',_display_command(command))
        return 0 if args.dry_run else subprocess.run(command,check=False).returncode
    except (PostprocessError,OSError,ValueError) as exc:
        print('[ERROR]',exc,file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
