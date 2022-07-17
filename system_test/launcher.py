'''
runs all tests under the given directory in parallel
'''
import tempfile, shutil, sys, time
'''
can't actually use this module directly since workload is determined by the number of MPI ranks
no the number of python processes
'''
import multiprocessing
from pathlib import Path
from subprocess import Popen, PIPE
import argparse


parser = argparse.ArgumentParser(description='Run M7 system tests in parallel')
parser.add_argument('mode', type=str, help='test or redefine benchmarks')
parser.add_argument('m7_exe', type=str, help='path to M7 binary')
parser.add_argument('--mpirun', type=str, help='path to mpirun executable', default='mpirun')
parser.add_argument('--nslot', type=int, help='total number of slots to make available for M7 jobs',
        default=multiprocessing.cpu_count())
parser.add_argument('--update_period', type=float, help='interval between queue polls (secs)', default=1.0)

args = parser.parse_args()
assert args.mode in ('benchmark', 'test')

queue_dir = Path(tempfile.mkdtemp())
print(f'reading job queue from {queue_dir}')

running_jobs = {}

def nslot_avail():
    navail = args.nslot - sum(v for k, v in running_jobs.items())
    assert navail >= 0
    return navail

def update():
    for job in queue_dir.glob('*'):
        if '.' in job.name:
            try: 
                name, nslot = job.name.split('.')
                nslot = int(nslot)
            except ValueError: 
                # invalid filename
                continue
            if nslot <= nslot_avail():
                job.unlink()
                running_jobs[name] = nslot
        else:
            if job.name in running_jobs:
                del running_jobs[job.name]

sys_test_dir = Path().absolute()

tests = []
procs = []
for f in sorted(Path().rglob('test.py')):
    tests.append({'path': f.parent})
    cmd = f'env PYTHONPATH={sys_test_dir} python {f.name} '
    cmd += f'{args.mode} {args.m7_exe} --mpirun={args.mpirun} --queue_dir={queue_dir}'
    procs.append(Popen(cmd, shell=True, cwd=f.parent, stdout=PIPE, stderr=PIPE))

while any(proc.poll() is None for proc in procs):
    time.sleep(args.update_period)
    update()

for iproc, proc in enumerate(procs):
    out, err = proc.communicate()
    tests[iproc]['out'], tests[iproc]['err'] = str(out, 'utf8'), str(err, 'utf8')

def result(test):
    if 'FAILURE' in test['err']: return 'FAILED'
    elif 'SKIP' in test['err']: return 'SKIPPED'
    return 'PASSED'


for test in tests:
    print('{:<50} {}'.format(str(test['path']), result(test)))

shutil.rmtree(queue_dir)
