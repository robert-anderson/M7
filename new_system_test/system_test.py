from subprocess import Popen, PIPE
from pathlib import Path
from sys import argv
import os, shutil, h5py, time, argparse
import numpy as np

import resource_manager

# test functionality is split into two concepts:
# - Static: verify current test instance against some external asset
# - Regression: verify current test instance against a previous benchmark instance

# set up paths for all relevant directories
# working directory: where the resource_manager keeps track of the number of ranks in use
WRK_DIR = Path(os.getcwd()).resolve()
# test definition: where the test script importing this file is located
DEF_DIR = Path(argv[0]).parent.resolve()
# benchmark: contains the artefacts defining a passing run
BMK_DIR = DEF_DIR/'benchmark'
# run: the working directory for instances of the tested program
RUN_DIR = DEF_DIR/'tmp'
# assets: static assets shared by multiple test definitions are defined under this dir
AST_DIR = Path(__file__).parent.parent/'assets'

# if this test does not have a defined benchmark, regression tests are skipped
HAVE_BMK = BMK_DIR.exists()
# if there is already a tmp directory, remove it
shutil.rmtree(RUN_DIR, ignore_errors=True)
RUN_DIR.mkdir()

assert AST_DIR.exists()

parser = argparse.ArgumentParser(description='Run this M7 system test')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('mpirun', type=str, help='path to mpirun executable', default='mpirun')
args = parser.parse_args()

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

def run(config_fname='config.yaml', nrank=1, copy_deps=[], link_deps=[]):
    cmd = f'{args.mpirun} -n {nrank} {args.m7_exe} {config_fname}'
    # config is copied so that it is retained when copied to benchmark
    copy_deps.append(config_fname)
    for dep in copy_deps: bring(dep, 'copy')
    for dep in link_deps: bring(dep, 'link')
    with resource_manager.instance(nrank):
        print(cmd)
        out, err = shell(cmd, RUN_DIR)
        assert not len(err), 'error stream non-empty'
