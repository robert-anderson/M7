from subprocess import Popen, PIPE
from pathlib import Path
import sys, os, shutil, h5py, argparse
import numpy as np
import pickle as pkl

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

# loop through root paths, if path is not found under any of these, then assume path is absolute
def resolve(root_order, path):
    for root in root_order:
        tmp = root/path
        if tmp.exists(): return tmp.resolve()
    if path.exists(): return path.resolve()
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
    src = resolve([DEF_DIR, AST_DIR], src_path)

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
    # 2: test not applicable to this binary
    # 3: no ref but test contains compare_ checks
    sys.exit(2 + bool(no_ref))

def fail(static, msg):
    # exit codes:
    # 4: static failure
    # 5: comparative failure
    print(f'{"STATIC" if static else "COMPARATIVE"} failure: {msg}')
    sys.exit(4 + bool(static))

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
                except IndexError: icolumn_end = self.ncolumn()
                return np.arange(icolumn_start, icolumn_end)
        return None

    def stats_columns(self, field_name_hint):
        return self.data[:, self.field_column_range(field_name_hint)]

def run(config_fname='config.yaml', nrank=1, copy_deps=[], link_deps=[]):
    cmd = f'{args.mpirun} -n {nrank} {args.m7_exe} {config_fname}'
    # config is copied so that it is retained when copied to ref
    copy_deps.append(config_fname)
    for dep in copy_deps: bring(dep, 'copy')
    for dep in link_deps: bring(dep, 'link')
    with resource_manager.instance(nrank):
        out, err = shell(cmd, RUN_DIR)
        assert not len(err), f'error stream non-empty: {err}'

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
    run = StatsFile(RUN_DIR/fname).stats_columns(field_name_hint)
    ref = StatsFile(REF_DIR/fname).stats_columns(field_name_hint)
    if not np.allclose(run, ref): fail(False, f'stats field "{field_name_hint}"')

def compare_nw(fname='M7.stats'): compare_stats_field('WF L1 norm', fname)
def compare_ref_weight(fname='M7.stats'): compare_stats_field('Reference weight', fname)
def compare_shift(fname='M7.stats'): compare_stats_field('Diagonal shift', fname)
def compare_ninit(fname='M7.stats'): compare_stats_field('Initiator', fname)
def compare_nocc_mbf(fname='M7.stats'): compare_stats_field('Occupied MBFs', fname)

def compare_rdm_archives(fname='M7.rdm.h5'):
    if not DO_COMPS: return

    run = h5py.File(RUN_DIR/fname, 'r')
    ref = h5py.File(REF_DIR/fname, 'r')
    for section in ('archive', 'spinfree'):
        if not section in ref.keys(): continue
        b = ref[section]
        r = run[section]
        keys = tuple(map(str, b.keys()))
        if set(r.keys()) != set(b.keys()):
            # different ranks of RDM accumulated than in benchmark
            fail(False, f'RDM "{section}" groups contain different keys')
        for key in keys:
            if key=='norm': continue
            if not np.array_equal(r[key]['indices'], b[key]['indices']):
                fail(False, f'index array of RDM {key}')
            if not np.allclose(np.array(r[key]['values']), np.array(b[key]['values'])): 
                fail(False, f'value array of RDM {key}')

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
    stats = StatsFile(RUN_DIR/fname).stats_columns(field_name_hint)
    mean, err = block(stats[-opts.npoint:], opts.nblock)
    err *= opts.err_scale
    if not within_error(ref_value, mean, err): 
        fail(True, f'stats field "{field_name_hint}" has mean and error {mean:.5e} +/- {err:.3e}, but ref is {ref_value:.5e}')

def check_shift(ref_value, fname='M7.stats', opts=BlockOpts()):
    check_stats_field(ref_value, 'Diagonal shift', fname, opts)

def check_proje(ref_value, fname='M7.stats', opts=BlockOpts()):
    check_stats_field(ref_value, 'Reference-projected energy', fname, opts)

def check_rdm_energy(ref_value, fname='M7.mae.stats', rtol=1e-5, atol=1e-8):
    stats = StatsFile(RUN_DIR/fname).stats_columns('Energy estimate from RDMs').ravel()
    mean = stats[-1]
    if not np.isclose(ref_value, mean, rtol, atol):
        fail(True, f'RDM energy has mean {mean:.5e}, but ref is {ref_value:.5e}')

def load_spinfree_hdf5_rdm(group):
    inds = np.array(group['indices'])
    values = np.array(group['values'])
    extent = max(inds.flatten())+1
    nind = inds.shape[1]
    rdm = np.zeros((extent,)*nind)
    for i, row in enumerate(inds): rdm[tuple(row)] = values[i]
    return rdm

def check_spinfree_rdms(h5_path, pkl_path, keys, tol=1e-5):
    h5_file = h5py.File(resolve([RUN_DIR], h5_path), 'r')
    pkl_fname = resolve([AST_DIR], pkl_path)
    with open(pkl_fname, 'rb') as f: pkl_rdms = pkl.load(f)

    for key in keys:
        h5_rdm = load_spinfree_hdf5_rdm(h5_file[f'spinfree/{key}'])
        pkl_rdm = pkl_rdms[key]
        abs_err = np.abs(h5_rdm - pkl_rdm)
        max_indices = np.unravel_index(np.argmax(abs_err), abs_err.shape)
        max_diff = abs_err[max_indices]
        if max_diff > tol:
            run_val = h5_rdm[max_indices]
            pkl_val = pkl_rdm[max_indices]
            fail(True, f'RDM {key} element {max_indices} value {run_val:.5e} does not equal reference {pkl_val:.5e} within tol {tol:.1e}')
