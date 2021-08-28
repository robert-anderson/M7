from pathlib import Path
import os, tempfile, shutil, importlib, sys, m7_test

try:
    m7_test.M7_EXE_PATH = sys.argv[1]
except IndexError:
    raise SystemExit('specify a path to the M7 executable')

if not m7_test.verify_m7_exe(): 
    raise SystemExit('path specified was not detected to be a valid M7 executable')

try:
    m7_test.MPIRUN_PATH = sys.argv[2]
except IndexError:
    raise SystemExit('specify a path to the mpirun executable')

try:
    rootdir = sys.argv[3]
    assert os.path.exists(rootdir), 'specify a valid root directory to run tests for'
except IndexError:
    rootdir = '.'

paths = tuple(Path(rootdir).rglob('__init__.py'))

home = os.getcwd()
dirnames = []
results = []
reasons = []
for path in paths:
    try: shutil.rmtree('./tmp')
    except FileNotFoundError: pass
    dirnames.append(os.path.dirname(path))
    print(f'running test "{dirnames[-1]}..."')
    shutil.copytree(dirnames[-1], './tmp')
    os.chdir('./tmp')
    spec = importlib.util.spec_from_file_location('*', '__init__.py')
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
        results.append(True)
        reasons.append('')
    except Exception as ex: 
        results.append(False)
        reasons.append(str(ex))
    os.chdir(home)


for dirname, result, reason in zip(dirnames, results, reasons):
    print('{:<30} {} {}'.format(dirname, 'PASSED' if result else 'FAILED:', reason))
