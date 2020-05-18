from subprocess import Popen, PIPE
import sys, os, re

m7_root = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+'/..')

long_opt_regex = re.compile(r'\-\-[a-z_]+')

def bytes_to_string(b):
    if (sys.version_info > (3, 0)):
        return str(b, 'utf8')
    else:
        return str(b)

def dict_iter(d):
    if (sys.version_info > (3, 0)):
        return d.items()
    else:
        return d.iteritems()


def find_similar(word, words):
    def similar(candidate):
        lword = word.lower()
        lcand = candidate.lower()
        return lword in lcand or lcand in lword
    return filter(similar, words)

class Options:
    legal_opts = {}
    def __init__(self, exepath):
        # run the executable with just the help option to determine the configurable options
        stdout, stderr = Popen('{} -h'.format(exepath), shell=1, stdout=PIPE, stderr=PIPE).communicate()
        strout = bytes_to_string(stdout)
        for line in strout.split('\n'):
            try:
                optname = long_opt_regex.search(line).group()
            except AttributeError:
                continue
            optname = optname[2:]
            setattr(self, optname, None)
            self.legal_opts[optname] = None

    def to_command_string(self):
        # get options from key-value pairs in this __dict__
        for k in self.__dict__.keys(): assert k in self.legal_opts, \
            'Unrecognized option "{}", perhaps you meant: \n\t{}'.format(k, '\n\t'.join(find_similar(k, self.legal_opts.keys())))
        return ' '.join(['--{}={}'.format(k, v) for k, v in 
                dict_iter(dict(filter(lambda kv : kv[1] is not None, dict_iter(self.__dict__))))])
            

class Instance:
    def __init__(self, nthread=1, nmpirank=1, build_dir='build', fnameout='M7.out', fnameerr='M7.err', debug=False, opts=None):
        self.nthread = nthread
        self.nmpirank = nmpirank
        self.fnameout = fnameout
        self.fnameerr = fnameerr
        self.exepath = os.path.abspath(m7_root+'/{}/src/{}'.format(build_dir, 'debug' if debug else 'release'))
        assert os.path.exists(self.exepath)
        self.opts = Options(self.exepath) if opts is None else opts

    def run(self):
        with open(self.fnameout, 'w') as fout, open(self.fnameerr, 'w') as ferr:
            cmd = 'env OMP_NUM_THREADS={} mpirun -n {} {} {}'.format(self.nthread, self.nmpirank, self.exepath, self.opts.to_command_string())
            fout.write(cmd+'\n')
            Popen(cmd, shell=1, stdout=fout, stderr=ferr).wait()

