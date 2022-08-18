'''
runs all tests under the given directory in parallel
'''
import tempfile, shutil, sys, time
'''
can't actually use this module directly since workload is determined by the number of MPI ranks
not the number of python processes
'''
import multiprocessing
from pathlib import Path
from subprocess import Popen, PIPE
import argparse


parser = argparse.ArgumentParser(description='Run M7 system tests in parallel')
parser.add_argument('mode', type=str, help='test or redefine benchmarks')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('asset_dir', type=str, help='path to asset directory')
parser.add_argument('--root_dir', type=str, help='path to directory below which all tests should be performed',
        default='.')
parser.add_argument('--mpirun', type=str, help='path to mpirun executable', default='mpirun')
parser.add_argument('--nslot', type=int, help='total number of slots to make available for M7 jobs',
        default=multiprocessing.cpu_count())
parser.add_argument('--update_period', type=float, help='interval between queue polls (secs)', default=1.0)

args = parser.parse_args()
assert args.mode in ('benchmark', 'test')

'''
set a temporary directory to act as a job queue
'''
queue_dir = Path(tempfile.mkdtemp())
print(f'reading job queue from {queue_dir}')

'''
jobs currently running.
keys are the uuid the jobs generate for themselves
values are the number of slots allocated to the job
'''
running_jobs = {}

def nslot_avail():
    navail = args.nslot - sum(v for k, v in running_jobs.items())
    assert navail >= 0
    return navail

def update():
    '''
    queued jobs will touch a file in the queue directory named with the uuid (name)
    followed by a '.', then the number of slots requested
    '''
    for job in queue_dir.glob('*'):
        if '.' in job.name:
            # this is a request by a job for resources
            try: 
                name, nslot = job.name.split('.')
                nslot = int(nslot)
            except ValueError: 
                # invalid filename
                continue
            if nslot <= nslot_avail():
                # remove the request
                job.unlink()
                running_jobs[name] = nslot
        else:
            # this is indicating that the job is completed
            if job.name in running_jobs:
                del running_jobs[job.name]

sys_test_dir = Path().absolute()

tests = []
# Popen process instances
procs = []
'''
find all directories under root_dir which contain a test.py script
'''
for f in sorted(Path(args.root_dir).rglob('test.py')):
    tests.append({'path': f.parent})
    cmd = f'env PYTHONPATH={sys_test_dir} python {f.name} '
    cmd += f'{args.mode} {args.m7_exe} {args.asset_dir} --mpirun={args.mpirun} --queue_dir={queue_dir}'
    procs.append(Popen(cmd, shell=True, cwd=f.parent, stdout=PIPE, stderr=PIPE))

'''
poll() method of process returns None if it is not complete
keep polling all processes until there are not any returning None
'''
while any(proc.poll() is None for proc in procs):
    time.sleep(args.update_period)
    update()

'''
all processes are completed, call communicate method to get stdout and stderr outputs
'''
for iproc, proc in enumerate(procs):
    out, err = proc.communicate()
    test = tests[iproc]
    test['out'], test['err'], test['code'] = str(out, 'utf8'), str(err, 'utf8'), proc.returncode

'''
based on stderr and return code, decide the test outcome
'''
def result(test):
    if 'FAILURE' in test['err']: return test['err'].split('Exception: ')[1]
    elif 'SKIP' in test['err']: return 'SKIPPED'
    elif 'No benchmark defined' in test['err']: return 'BENCHMARK UNDEFINED'
    elif test['code']: return 'TEST ERROR:'+test['err']
    return 'NEW BENCHMARK DEFINED' if (args.mode=='benchmark') else 'PASSED'

for test in tests:
    print('{:<50} {}'.format(str(test['path']), result(test)))

shutil.rmtree(queue_dir)
