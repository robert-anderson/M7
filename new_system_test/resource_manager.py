import os, time
from contextlib import contextmanager

NRANK_AVAIL = os.cpu_count()
POLL_PERIOD_SECS = 0.25
INSTANCE_COUNT_FNAME = '.instance_count'

def read_ninstance():
    try:
        with open(INSTANCE_COUNT_FNAME, 'r') as f: return int(f.read())
    except (ValueError, FileNotFoundError):
        return 0

def update_ninstance(n):
    with open(INSTANCE_COUNT_FNAME, 'w') as f: f.write(str(int(n)))

def poll_until(fn, *args):
    while True:
        if fn(*args): break
        time.sleep(POLL_PERIOD_SECS)

def acquire_lock_file():
    try:
        with open('.lock', "x") as fn: 
            return True
    except FileExistsError:
        pass
    return False

@contextmanager
def lock():
    acquired = False
    try:
        poll_until(acquire_lock_file)
        acquired = True
        yield
    finally:
        if (acquired): os.remove('.lock')

def acquire_instance(nrank):
    with lock():
        try: ncurrent = read_ninstance()
        except FileNotFoundError: ncurrent = 0
        if ncurrent + nrank <= NRANK_AVAIL:
            update_ninstance(ncurrent + nrank)
            return True
    return False

@contextmanager
def instance(nrank, on_exit=None):
    assert nrank <= NRANK_AVAIL, 'insufficient total resources for this job'
    acquired = False
    try:
        poll_until(acquire_instance, nrank)
        acquired = True
        yield
    finally:
        if (acquired):
            with lock():
                ncurrent = read_ninstance()
                update_ninstance(ncurrent - nrank)
            if on_exit is not None: on_exit()

if __name__=='__main__':
    # test with a fake task (sleep for nsec seconds) on nrank ranns
    from sys import argv
    nsec = int(argv[1])
    nrank = int(argv[2])

    with instance(nrank):
        time.sleep(nsec)
