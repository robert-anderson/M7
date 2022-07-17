from subprocess import Popen, PIPE
from pathlib import Path
from sys import argv
import tempfile, os, shutil
import numpy as np
from uuid import uuid4 as uuid
import time, argparse


parser = argparse.ArgumentParser(description='Run this M7 system test')
parser.add_argument('mode', type=str, help='test or redefine benchmarks')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('--mpirun', type=str, help='path to mpirun executable', default=None)
parser.add_argument('--queue_dir', type=str, help='path queue for parallel testing', default=None)
args = parser.parse_args()

assert args.mode in ('benchmark', 'test')

def_dir = Path()
run_dir = def_dir/'tmp'
shutil.rmtree(run_dir, ignore_errors=True)
run_dir.mkdir()

bench_dir = def_dir/'benchmark'

compile_defs = None
benchmarking = args.mode=='benchmark'

# if testing, make sure we have something to compare against
if not benchmarking: assert bench_dir.exists()

request_delay = 0.2

def shell(cmd, cwd=Path()):
    tmp = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, cwd=cwd).communicate()
    return str(tmp[0], 'utf8'), str(tmp[1], 'utf8')

for line in shell(args.m7_exe)[0].split('\n'):
    if 'compile definitions' in line: compile_defs = {}
    elif 'Input specification' in line: break
    if compile_defs is None: continue
    split = line.strip().split('|')
    if len(split)==4: compile_defs[split[1].strip()] = split[2].strip()

# validate given paths
assert compile_defs is not None, "invalid M7 binary"
if args.mpirun is not None: assert 'mpi' in shell(f'{args.mpirun} --version')[0].lower(), 'invalid mpirun path'

class StatsFile:
    columns = []
    data = None
    def __init__(self, fname):
        with open(fname, 'r') as f:
            for line in f.readlines():
                if not line.startswith('#'): break
                split = line[1:].strip().split('.')
                try: i = int(split[0])-1
                except ValueError: continue
                self.columns.append((i, split[1].split('(')[0].strip()))
        self.data = np.loadtxt(fname)
    def lookup_column(self, start):
        for column in self.columns:
            if column[1].lower().startswith(start.lower()): return column
        return None

class TestInstance:
    error_msg = ""
    irun = 0
    irun_last_stats_update = 0
    run_stats = None
    bench_stats = None

    def stats(self, fname='M7.stats'):
        if benchmarking or len(self.error_msg): return
        if self.irun_last_stats_update != self.irun:
            self.run_stats = StatsFile(run_dir/fname)
            self.bench_stats = StatsFile(bench_dir/fname)
            self.irun_last_stats_update = self.irun
            if self.run_stats.columns != self.bench_stats.columns:
                self.error_msg = "incompatible column names"
        return self.run_stats, self.bench_stats

    def __del__(self):
        print(self.error_msg)
        if not benchmarking: return
        shutil.rmtree(bench_dir, ignore_errors=True)
        shutil.move(run_dir, bench_dir)


# when this goes out of scope, the benchmark will be adjusted if requested
instance = TestInstance()
def fail(msg):
    raise Exception(f'FAILURE: {msg}')

def skip():
    raise Exception('SKIP')

def stats_columns(col_name, fname='M7.stats'):
    stats = instance.stats(fname)
    column = stats[0].lookup_column(col_name)
    assert column is not None
    return stats[0].data[:,column[0]], stats[1].data[:,column[0]]

def bring(fname, src_dir=None):
    src = Path()/fname if src_dir is None else Path(src_dir)/fname
    dst = run_dir/fname
    assert src.exists()
    assert not dst.exists()
    os.symlink(src.absolute(), dst.absolute())

def run(config='config.yaml', nrank=1, assets=tuple()):
    assert nrank
    if not (def_dir/config).exists(): config = ''
    for asset in assets+(config,): bring(asset)
    cmd = f'{args.m7_exe} {config}' 
    if args.mpirun is not None: cmd = f'{args.mpirun} -n {nrank} '+cmd
    else: assert nrank==1, 'require more than one rank but no mpirun path is given'

    if args.queue_dir is not None:
        job_name = uuid()
        request_path = Path(args.queue_dir)/f'{job_name}.{nrank}'
        done_path = Path(args.queue_dir)/f'{job_name}'
        request_path.touch()
        while request_path.exists(): time.sleep(request_delay)

    instance.irun+=1
    shell(cmd, run_dir)

    if args.queue_dir is not None: done_path.touch()

    if benchmarking:
        for asset in assets+(config,): os.unlink(run_dir/asset)


def skip_if(compile_def, value):
    assert compile_def in compile_defs
    if compile_defs[compile_def]==value: skip()

def skip_unless(compile_def, value):
    assert compile_def in compile_defs
    if compile_defs[compile_def]!=value: skip()

def check_nw(fname='M7.stats'):
    if benchmarking: return
    run, bench = stats_columns('WF L1 norm', fname)
    if not np.allclose(run, bench): fail('walker trajectories do not agree')
