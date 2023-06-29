from subprocess import Popen
from pathlib import Path
import sys
import os, shutil, h5py, time, argparse
import numpy as np

# this script simply dispatches a number of test scripts

parser = argparse.ArgumentParser(description='Run a suite of M7 system tests')
parser.add_argument('mode', type=str, choices=['test', 'redef'],
        help='"test" to run tests or "redef" to run tests and redefine the refs of the static-passing ones')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('mpirun', type=str, help='path to mpirun executable', default='mpirun', nargs='?')
parser.add_argument('-p', '--paths', nargs='+', default=[])
args = parser.parse_args()

print(f'Total MPI slots: {os.cpu_count()}')

for path in args.paths: assert Path(path).exists(), path

procs = []
for path in args.paths:
    cmd = [sys.executable, Path(path).resolve(), args.m7_exe, args.mpirun, str(int(args.mode=='redef'))]
    procs.append(Popen(cmd, cwd=os.getcwd()))

# None = not finished
failed = [None for p in procs]

def colored(text, name):
    try:
        code = dict(
            BLACK = '\033[30m',
            RED = '\033[31m',
            GREEN = '\033[32m',
            YELLOW = '\033[33m',
            BLUE = '\033[34m',
            MAGENTA = '\033[35m',
            CYAN = '\033[36m',
            WHITE = '\033[37m'
        )[name.upper()]
    except KeyError: return name
    return f'{code}{text}\033[0m'

def outcome_string(exit_code, path):
    path = str(path)
    if exit_code < 1: return colored('PASS: '+path, 'green')
    elif exit_code < 2: return colored('SCRIPT ERROR: '+path, 'red')
    elif exit_code < 4: return colored('SKIP: '+path, 'yellow')
    else: return colored('FAIL: '+path, 'cyan')

from resource_manager import poll_until, read_ninstance
def all_procs_done():
    ndone = 0
    for i, proc in enumerate(procs):
        if failed[i] is not None: 
            ndone += 1
            continue
        failed[i] = proc.poll()
        if failed[i] is not None:
            ndone += 1
            print(outcome_string(failed[i], Path(args.paths[i]).parent))
            print(f'MPI slots in use: {read_ninstance()}/{os.cpu_count()}')
    return ndone == len(failed)

poll_until(all_procs_done)

def redef(script_path):
    def_dir = Path(script_path).resolve().parent
    run_dir = def_dir / 'tmp'
    ref_dir = def_dir / 'ref'
    # delete all symlinks in run directory
    for filename in os.listdir(run_dir):
        file_path = os.path.join(run_dir, filename)
        if os.path.islink(file_path): os.remove(file_path)
    shutil.rmtree(ref_dir, ignore_errors=True)
    shutil.move(run_dir, ref_dir)

assert not any(failed), 'not all tests passed'
if (args.mode=='redef'):
    print('redefining references for all statically-passing tests:')
    for i, exit_code in enumerate(failed):
        if exit_code is None: continue
        print(args.paths[i])
        redef(args.paths[i])
        
