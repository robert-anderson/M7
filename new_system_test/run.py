from subprocess import Popen, PIPE
from pathlib import Path
import sys
import os, shutil, h5py, time, argparse
import numpy as np

# this script simply dispatches a number of test scripts

parser = argparse.ArgumentParser(description='Run a suite of M7 system tests')
#parser.add_argument('--rebench', type=bool, help='path to M7 binary')
#parser.add_argument('--feature', action='store_true')
#parser.add_argument('--no-feature', dest='feature', action='store_false')
#parser.set_defaults(feature=True)

parser.add_argument('mode', type=str, choices=['test', 'rebench'],
        help='"test" to run tests or "rebench" to run tests and rebenchmark the static-passing ones')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('mpirun', type=str, help='path to mpirun executable', default='mpirun')
parser.add_argument('-p', '--paths', nargs='+', default=[])
args = parser.parse_args()

for path in args.paths: assert Path(path).exists()

procs = []
for path in args.paths:
    cmd = [sys.executable, Path(path).resolve(), args.m7_exe, args.mpirun]
    procs.append(Popen(cmd, cwd=os.getcwd()))

# None = not finished
failed = [None for p in procs]

def all_procs_done():
    ndone = 0
    for i, proc in enumerate(procs):
        if failed[i] is not None:
            ndone += 1
            print(('PASSED' if not failed[i] else 'FAILED')+f': {args.paths[i]}')
            continue
        failed[i] = proc.poll()
    return ndone == len(failed)

from resource_manager import poll_until
poll_until(all_procs_done)

def rebench(script_path):
    def_dir = Path(script_path).resolve().parent
    run_dir = def_dir / 'tmp'
    bmk_dir = def_dir / 'benchmark'
    # delete all symlinks in run directory
    for filename in os.listdir(run_dir):
        file_path = os.path.join(run_dir, filename)
        if os.path.islink(file_path): os.remove(file_path)
    shutil.rmtree(bmk_dir, ignore_errors=True)
    shutil.move(run_dir, bmk_dir)

if (args.mode=='rebench'):
    print('rebenchmarking all statically-passing tests:')
    for i, exit_code in enumerate(failed):
        if exit_code: continue
        print(args.paths[i])
        rebench(args.paths[i])
        
assert not any(failed), 'not all tests passed'
