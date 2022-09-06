import itertools
import math, sys
import numpy as np

'''
the spin signature (spinsig) is given by a bitstring created from the spins of each string (creation / annihilation)
the spincases are a series of contiguous labels given to the valid (spin conserving) pairs of spinsigs
'''

def to_file(name_data_signed, header_lines, fname, guard_var=None):
    header_lines.insert(0, '')
    header_lines.insert(0, 'this file was generated programmatically and should not be modified by hand.')
    header_lines = ['\n * '+line for line in header_lines]
    if guard_var is None: guard_var = 'M7_'+fname.split('.')[0].upper()+'_'+fname.split('.')[-1].upper()

    basic_fmt_fn = lambda i: str(i)

    with open(fname, 'w') as f:
        f.write('/*'+(''.join(header_lines))+'\n */\n')
        f.write(f'#ifndef {guard_var}\n')
        f.write(f'#define {guard_var}\n')
        f.write('#include "M7_lib/defs.h"\n')
        for name, data, signed in name_data_signed:
            fmt_fn = basic_fmt_fn
            if not signed: fmt_fn = lambda i: basic_fmt_fn(i) if i>=0 else '~0ul'
            np.set_printoptions(threshold=sys.maxsize, formatter={'int': fmt_fn}, linewidth=120)
            content = data.__repr__()[6:-1].replace('[', '{').replace(']', '}')
            dims_str = ''.join(f'[{i}]' for i in data.shape)
            f.write(f'constexpr {"uint_t" if not signed else "int"} {name}{dims_str} = \n      {content};\n\n')
        f.write(f'#endif //{guard_var}\n')

def n_set_bits(n):
    if n<1: return 0
    tot = 0
    for i in range(int(math.log(n, 2.0))+1):
        if (n>>i)&1: tot+=1
    return tot

def int_to_tup(n, rank):
    return tuple((n>>i)&1 for i in range(rank))

class SpinSigPairs(list):
    def __init__(self, maxrank):
        self.maxrank = maxrank
        list.__init__(self)
        for i in range(2**self.maxrank):
            for j in range(i):
                if n_set_bits(i)-n_set_bits(j): continue
                self.append([i, j])
                self.append([j, i])
            self.append([i, i])

if __name__=='__main__':
    maxrank = 4
    arr = np.zeros((2**maxrank,)*2, dtype=int)
    arr[:,:] = -1
    for icase, (ihole, ielec) in enumerate(SpinSigPairs(maxrank)):
        arr[ihole, ielec] = icase

    
    header_lines = [
        'only a subset of spin signature pairs are valid for RDM processing, so the',
        'array in this file maps the creation and annihilation spin signature to a spin',
        'signature "case_id" for subsequent processing']

    to_file([('spinsig_map', arr, False)], header_lines, 'SpinsigMap.h')


