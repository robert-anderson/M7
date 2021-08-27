from pathlib import Path
import os, tempfile, shutil, importlib

paths = tuple(Path('.').rglob('__init__.py'))

home = os.getcwd()
dirnames = []
results = []
reasons = []
for path in paths:
    dirnames.append(os.path.dirname(path))
    shutil.copytree(dirnames[-1], 'tmp')
    shutil.move('tmp', dirnames[-1])
    os.chdir(dirnames[-1]+'/tmp')
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
    shutil.rmtree(dirnames[-1]+'/tmp')


for dirname, result, reason in zip(dirnames, results, reasons):
    print('{:<30} {} {}'.format(dirname, 'PASSED' if result else 'FAILED:', reason))
