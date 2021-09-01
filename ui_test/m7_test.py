'''
just two module-level constants defining the path to the mpirun command are needed, the rest of this
file contains functions to invoke M7 instances and verify their output
'''
from subprocess import Popen, PIPE
import pandas as pd
import os

M7_EXE_PATH = None
MPIRUN_PATH = None


def run_command(cmd):
    stdout, stderr = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
    return str(stdout, 'utf8'), str(stderr, 'utf8')

def verify_m7_exe():
    if M7_EXE_PATH is None: return False
    stdout, stderr = run_command(M7_EXE_PATH)
    return 'Configuration document "FCIQMC options"' in stdout.split('\n')[1]

def run(yaml_path='config.yaml', nrank=1):
    assert os.path.exists(yaml_path), 'yaml path is invalid'
    assert M7_EXE_PATH is not None and MPIRUN_PATH is not None, \
            'invalid paths, please run tests through the driver script'
    return run_command('{} -n {} {} {}'.format(MPIRUN_PATH, nrank, M7_EXE_PATH, yaml_path))
