from subprocess import Popen, PIPE

M7_exe = "/home/rja/CLionProjects/M7/cmake-build-debug/bin/M7"

stdout, stderr = Popen(M7_exe, stdout=PIPE, stderr=PIPE, shell=True).communicate()

lines = str(stdout, 'utf8').split('\n')


def get_compile_defs(lines):
    d = None
    for line in lines:
        if 'compile definitions' in line: d = {}
        elif 'Input specification' in line: return d
        if d is None: continue
        split = line.strip().split('|')
        if len(split)==4: d[split[1].strip()] = split[2].strip()

print(get_compile_defs(lines))
