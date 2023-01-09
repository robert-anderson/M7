from subprocess import Popen, PIPE
from pathlib import Path
from sys import argv
import tempfile, os, shutil, h5py, time, argparse
import numpy as np
from uuid import uuid4 as uuid


parser = argparse.ArgumentParser(description='Run this M7 system test')
parser.add_argument('mode', type=str, help='test or redefine benchmarks')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('asset_dir', type=str, help='path to asset directory')
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
if not benchmarking: assert bench_dir.exists(), 'No benchmark defined'

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

asset_dir = Path(args.asset_dir)
assert asset_dir.exists()

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

    irun_last_rdms_update = 0
    run_rdms = None
    bench_rdms = None


    redefine_benchmark = True

    def stats(self, fname='M7.stats'):
        if len(self.error_msg): return
        if self.irun_last_stats_update != self.irun:
            self.irun_last_stats_update = self.irun
            self.run_stats = StatsFile(run_dir/fname)
            if not benchmarking: 
                self.bench_stats = StatsFile(bench_dir/fname)
                if self.run_stats.columns != self.bench_stats.columns:
                    self.error_msg = "incompatible column names"
        return self.run_stats, self.bench_stats

    def rdms(self, fname='M7.h5'):
        if len(self.error_msg): return
        if self.irun_last_rdms_update != self.irun:
            self.irun_last_rdms_update = self.irun
            self.run_rdms = h5py.File(run_dir/fname, 'r')
            if not benchmarking: 
                self.bench_rdms = h5py.File(bench_dir/fname, 'r')
        return self.run_rdms, self.bench_rdms

    def __del__(self):
        print(self.error_msg)
        if not benchmarking or not self.redefine_benchmark: return
        shutil.rmtree(bench_dir, ignore_errors=True)
        shutil.move(run_dir, bench_dir)


# when this goes out of scope, the benchmark will be adjusted if requested
instance = TestInstance()
def fail(msg):
    instance.redefine_benchmark = False
    raise Exception(f'FAILURE: {msg}')

def skip():
    instance.redefine_benchmark = False
    raise Exception('SKIP')

def stats_columns(col_name, fname='M7.stats'):
    stats = instance.stats(fname)
    column = stats[0].lookup_column(col_name)
    assert column is not None
    if not benchmarking:
        return stats[0].data[:,column[0]], stats[1].data[:,column[0]]
    else:
        return stats[0].data[:,column[0]], None

def is_vector(obj):
    return isinstance(obj, tuple) or isinstance(obj, list)

def resolve_src_path(asset):
    if is_vector(asset):
        path = asset[0]
    else:
        path = asset
    # first, attempt to resolve path locally
    src = Path()/path
    # if this fails, look relative to the asset directory
    if not src.exists(): src = asset_dir.joinpath(path)
    assert src.exists(), f'path "{src}" found neither locally nor under assets directory'
    return src

def make_local_name(asset):
    # asset can be tuple
    if is_vector(asset):
        # return given local name
        name = asset[1]
    else:
        # inherit local name from that of src
        name = Path(asset).name
    return Path(run_dir/name)

def bring(asset):
    src = resolve_src_path(asset)
    dst = make_local_name(asset)
    assert not dst.exists()
    os.symlink(src.absolute(), dst.absolute())

def run(config='config.yaml', nrank=1, assets=[]):
    assert nrank
    if not (def_dir/config).exists(): config = ''
    else: shutil.copy(def_dir/config, run_dir/config)

    for asset in assets: bring(asset)
    cmd = f'{args.m7_exe} {config}' 
    if args.mpirun is not None: cmd = f'{args.mpirun} -n {nrank} '+cmd
    else: assert nrank==1, 'require more than one rank but no mpirun path is given'

    if args.queue_dir is not None:
        # we're using a queue
        # generate a unique name for the job
        job_name = uuid()
        request_path = Path(args.queue_dir)/f'{job_name}.{nrank}'
        # when this path comes into existence, the job is completed
        done_path = Path(args.queue_dir)/f'{job_name}'
        # request resources from the queue by creating a blank file at the request_path
        request_path.touch()
        # the request must wait its turn
        while request_path.exists(): 
            # deletion of the request indicates that the queue approves the job to run
            time.sleep(request_delay)

    instance.irun+=1
    out, err = shell(cmd, run_dir)
    if len(err): instance.redefine_benchmark = False

    # let the queue know that this process is done
    if args.queue_dir is not None: done_path.touch()

    if benchmarking and instance.redefine_benchmark:
        '''
        keep benchmarks clean of symlinks, since these are in general directory structure-dependent
        only unlink assets if this in benchmarking mode, otherwise leave links in place for debugging
        '''
        for asset in assets: os.unlink(make_local_name(asset))

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

def check_ninit(fname='M7.stats'):
    if benchmarking: return
    run, bench = stats_columns('Initiator', fname)
    if not np.allclose(run, bench): fail('initiator number trajectories do not agree')

def check_rdm_archives(fname='M7.rdm.h5'):
    if benchmarking: return
    run, bench = instance.rdms(fname)
    for section in ('archive', 'spinfree'):
        if not section in bench.keys(): continue
        b = bench[section]
        r = run[section]
        keys = tuple(map(str, b.keys()))
        assert set(r.keys())==set(b.keys()), 'different ranks of RDM accumulated than in benchmark'
        for key in keys:
            if key=='norm': continue
            if not np.array_equal(r[key]['indices'], b[key]['indices']):
                fail(f'index array of RDM {key} does not agree with benchmark')
            if not np.allclose(np.array(r[key]['values']), np.array(b[key]['values'])): 
                fail(f'value array of RDM {key} does not agree with benchmark')

'''
perform crude removal of serial correlation
'''
def block(series, nblock):
    avgs = np.zeros(nblock)
    blocklen = len(series)//nblock
    for i in range(nblock): avgs[i] = np.mean(series[i*blocklen:(i+1)*blocklen])
    return np.mean(series), np.sqrt(np.var(avgs)/nblock)

def within_error(correct_value, mean, error):
    return (correct_value > mean - error) and (correct_value <= mean + error)

'''
check that a stats column is statistically correct (within errorbars)
'''
def check_stat_correct(correct_value, column_name, nburnin, nblock=32, fname='M7.stats'):
    stats, _ = stats_columns(column_name, fname)
    mean, err = block(stats[nburnin:], nblock)
    if not within_error(correct_value, mean, err):
        fail(f'mean of column "{column_name}" ({mean:.6e}) is not correct ({correct_value:.6e}) within error bars (+/-{err:.3e})')

def check_shift_correct(correct_value, nburnin, nblock=32, fname='M7.stats'):
    check_stat_correct(correct_value, 'Diagonal shift', nburnin, nblock, fname)

def check_proje_correct(correct_value, nburnin, nblock=32, fname='M7.stats'):
    check_stat_correct(correct_value, 'Reference-projected energy', nburnin, nblock, fname)
