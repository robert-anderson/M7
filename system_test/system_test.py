from subprocess import Popen, PIPE
from pathlib import Path
import sys, os, shutil, h5py, argparse
import numpy as np

import resource_manager

# test functionality is split into two concepts:
# - Static: verify current test instance against some external asset
# - Comparative: verify current test instance against a trusted reference instance

# set up paths for all relevant directories
# working directory: where the resource_manager keeps track of the number of ranks in use
WRK_DIR = Path(os.getcwd()).resolve()
# test definition: where the test script importing this file is located
DEF_DIR = Path(sys.argv[0]).parent.resolve()
# reference: contains the artefacts defining a passing run
REF_DIR = DEF_DIR/'ref'
# run: the working directory for instances of the tested program
RUN_DIR = DEF_DIR/'tmp'
# assets: static assets shared by multiple test definitions are defined under this dir
AST_DIR = Path(__file__).parent.parent/'assets'

# if there is already a tmp directory, remove it
shutil.rmtree(RUN_DIR, ignore_errors=True)
RUN_DIR.mkdir()

assert AST_DIR.exists()

parser = argparse.ArgumentParser(description='Run this M7 system test')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('mpirun', type=str, help='path to mpirun executable', default='mpirun', nargs='?')
parser.add_argument('static_only', type=int, help='if non-zero, any defined reference is ignored', default=0, nargs='?')
args = parser.parse_args()

# if this test does not have a defined reference or static only is specified, comparative tests are skipped
DO_COMPS = REF_DIR.exists() and not bool(args.static_only)

# root/path is given priority. if it doesn't exist then assume path is absolute
def resolve_path(root, path):
    tmp = root/path
    if tmp.exists(): return tmp.resolve()
    return None

def is_vector(obj):
    return isinstance(obj, tuple) or isinstance(obj, list)

def make_local_name(asset):
    # asset can be tuple
    if is_vector(asset):
        # return given local name
        name = asset[1]
    else:
        # inherit local name from that of src
        name = Path(asset).name
    return Path(run_dir/name)

# "bring" a file on which the test depends into the RUN_DIR either by copy or soft symlink
def bring(path_or_pair, kind):
    assert kind in ('copy', 'link'), 'invalid dependency kind'
    try:
        src_path, dst_path = path_or_pair
    except (ValueError, TypeError):
        src_path, dst_path = path_or_pair, Path(path_or_pair).name

    # dst_path is always relative to the temporary run directory
    dst_path = RUN_DIR/dst_path
    # first look in this test's definition directory
    src = resolve_path(DEF_DIR, src_path)
    # if not found, try the assets directory
    if src is None:
        src = resolve_path(AST_DIR, src_path)
    assert src is not None, f'file dependency "{src_path}" not found'
    dst = Path(dst_path).resolve()
    if dst.exists(): os.unlink(dst)
    if kind=='copy':
        shutil.copy(src, dst)
    elif kind=='link':
        os.symlink(src, dst)

def shell(cmd, wd):
    tmp = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, cwd=wd).communicate()
    return str(tmp[0], 'utf8'), str(tmp[1], 'utf8')

COMPILE_DEFS = None
for line in shell(args.m7_exe, '.')[0].split('\n'):
    if 'compile definitions' in line: COMPILE_DEFS = {}
    elif 'Input specification' in line: break
    if COMPILE_DEFS is None: continue
    split = line.strip().split('|')
    if len(split)==4: COMPILE_DEFS[split[1].strip()] = split[2].strip()

# validate given paths
assert COMPILE_DEFS is not None, "invalid M7 binary"
MBF_TYPE = COMPILE_DEFS['many-body basis function'].split()[0]
MBF_TYPES = ('fermion', 'boson', 'fermion-boson')
HAM_ARITH = COMPILE_DEFS['Hamiltonian arithmetic'].split()[0]
HAM_ARITHS = ('real', 'complex')

def skip(no_ref):
    # exit codes:
    # 1: test not applicable to this binary
    # 2: no ref but test contains compare_ checks
    sys.exit(1 + bool(no_ref))

def fail(static):
    # exit codes:
    # 3: static failure
    # 4: comparative failure
    sys.exit(3 + bool(static))

def require_mbf_type(s):
    assert s.lower() in MBF_TYPES
    if s.lower()!=MBF_TYPE: skip(False)

def require_ham_arith(s):
    assert s.lower() in HAM_ARITHS
    if s.lower()!=HAM_ARITH: skip(False)

class StatsFile:
    fields = []
    data = None
    def __init__(self, fname):
        with open(fname, 'r') as f:
            for line in f.readlines():
                if not line.startswith('#'): break
                split = line[1:].strip().split('.')
                try: i = int(split[0])-1
                except ValueError: continue
                self.fields.append((i, split[1].split('(')[0].strip()))
        self.data = np.loadtxt(fname)

    def ncolumn(self):
        return self.data.shape[1]

    def field_column_range(self, field_name_hint):
        for i, field in enumerate(self.fields):
            if field[1].lower().startswith(field_name_hint.lower()): 
                icolumn_start = field[0]
                try: icolumn_end = self.fields[i+1][0]
                except IndexError: icolumn_end = ncolumn()
                return np.arange(icolumn_start, icolumn_end)
        return None

    def stats_columns(self, field_name_hint):
        return self.data[:, self.field_column_range(field_name_hint)]

ref_stats_file = StatsFile(REF_DIR/'M7.stats') if DO_COMPS else None
run_stats_file = None

def run(config_fname='config.yaml', nrank=1, copy_deps=[], link_deps=[]):
    cmd = f'{args.mpirun} -n {nrank} {args.m7_exe} {config_fname}'
    # config is copied so that it is retained when copied to ref
    copy_deps.append(config_fname)
    for dep in copy_deps: bring(dep, 'copy')
    for dep in link_deps: bring(dep, 'link')
    with resource_manager.instance(nrank):
        out, err = shell(cmd, RUN_DIR)
        assert not len(err), 'error stream non-empty'

    # update stats to those of this run
    global run_stats_file
    run_stats_file = StatsFile(RUN_DIR/'M7.stats')

def stats_columns(col_name, fname='M7.stats'):
    stats = instance.stats(fname)
    column = stats[0].lookup_column(col_name)
    assert column is not None
    if not benchmarking:
        return stats[0].data[:,column[0]], stats[1].data[:,column[0]]
    else:
        return stats[0].data[:,column[0]], None

# compare_ methods involve verification against the contents of the ref directory
def compare_stats_field(field_name_hint, fname='M7.stats'):
    if not DO_COMPS: return
    run = run_stats_file.stats_columns(field_name_hint)
    ref = ref_stats_file.stats_columns(field_name_hint)
    if not np.allclose(run, ref): fail(False)

def compare_nw(fname='M7.stats'): compare_stats_field('WF L1 norm', fname)
def compare_ref_weight(fname='M7.stats'): compare_stats_field('Reference weight', fname)
def compare_shift(fname='M7.stats'): compare_stats_field('Diagonal shift', fname)
def compare_ninit(fname='M7.stats'): compare_stats_field('Initiator', fname)
def compare_nocc_mbf(fname='M7.stats'): compare_stats_field('Occupied MBFs', fname)


'''
perform crude removal of serial correlation
'''
def block(series, nblock):
    avgs = np.zeros(nblock)
    blocklen = len(series)//nblock
    for i in range(nblock): avgs[i] = np.mean(series[i*blocklen:(i+1)*blocklen])
    return np.mean(series), np.sqrt(np.var(avgs)/nblock)

def within_error(ref_value, mean, error):
    return (ref_value > mean - error) and (ref_value <= mean + error)

class BlockOpts:
    def __init__(self, npoint=1000, nblock=32, err_scale=2.0):
        self.npoint, self.nblock, self.err_scale = npoint, nblock, err_scale

'''
check that a stats column is statistically correct (within errorbars)
'''
def check_stats_field(ref_value, field_name_hint, fname='M7.stats', opts=BlockOpts()):
    stats = run_stats_file.stats_columns(field_name_hint)
    mean, err = block(stats[-opts.npoint:], opts.nblock)
    err *= opts.err_scale
    if not within_error(ref_value, mean, err): fail(True)

def check_shift(ref_value, fname='M7.stats', opts=BlockOpts()):
    check_stats_field(ref_value, 'Diagonal shift', fname, opts)

def check_proje(ref_value, fname='M7.stats', opts=BlockOpts()):
    check_stats_field(ref_value, 'Reference-projected energy', fname, opts)


