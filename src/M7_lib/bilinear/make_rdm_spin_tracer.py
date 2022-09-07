import itertools
import math, sys
import numpy as np

'''
the spin signature (spinsig) is given by a bitstring created from the spins of each string (creation / annihilation)
the spincases are a series of contiguous labels given to the valid (spin conserving) pairs of spinsigs
'''


'''
an item shall be a dict with the following keys:
    name, data, signed, comment, dimnames
'''

def to_file(namespace, items, header_lines, fname, guard_var=None):
    consts = {}
    for item in items:
        assert len(item['dimnames'])==len(item['data'].shape)
        for symbol, val in zip(item['dimnames'], item['data'].shape):
            if symbol not in consts: consts[symbol] = val
            else: assert val==consts[symbol], f'incompatible values given for "{symbol}"'

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

        f.write(f'namespace {namespace} {{\n')
        for symbol, val in consts.items():
            f.write(f'    constexpr uint_t {symbol} = {val};\n')

        f.write('\n')
        for item in items:
            name = item['name']
            data = item['data']
            signed = item['signed']
            comment = item['comment']
            dimnames = item['dimnames']
            fmt_fn = basic_fmt_fn
            if not signed: fmt_fn = lambda i: basic_fmt_fn(i) if i>=0 else '~0ul'
            np.set_printoptions(threshold=sys.maxsize, formatter={'int': fmt_fn}, linewidth=120)
            content = data.__repr__()[6:-1].replace('[', '{').replace(']', '}')
            dims_str = ''.join(f'[{dimname}]' for dimname in dimnames)
            f.write(f'    constexpr {"uint_t" if not signed else "int"} {name}{dims_str} = \n      {content};\n\n')
        f.write(f'\n}} // namespace {namespace}\n')
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

items = []
for rank in (2, 3):
    arr = np.zeros((2**rank,)*2, dtype=int)
    arr[:,:] = -1
    for icase, (ihole, ielec) in enumerate(SpinSigPairs(rank)):
        arr[ihole, ielec] = icase


    #header_lines = [
    #    'only a subset of spin signature pairs are valid for RDM processing, so the',
    #    'array in this file maps the creation and annihilation spin signature to a spin',
    #    'signature "case_id" for subsequent processing']

    items.append({
        'name': f'c_spinsig_pairs_{rank}', 
        'data': arr, 
        'signed': False, 
        'comment': f'maps the creation and annihilation spinsigs to a pair-spinsig id',
        'dimnames': [f'c_nspinsig_{rank}',]*2})

'''
a "case" is identified by:
    - a spin signature case
    - a spatial signature case for holes
    - a spatial signature case for elecs

the hole and elec spin signatures are combined into a single case because
spin symmetry leads to a sparse matrix of valid combinations. see make_spinsig_map.py

the spatial signatures in principle have the same sparse property but point group
symmetry is only realised at runtime, so this is a memory optimisation which can't
be made. The cost in memory is minimal and the CPU cost is zero.

the spatial signature is an integer describing the spatial structure of an RDM element
if the RDM element has rank r, the spatial signature will be encoded in r-1 bits

example: if the elec part of the drawn spinned element is 1, 3, 4, 7
the spatial orbitals are                                  1, 2, 2, 4
and the bit representation would be                          0, 1, 0 (= 2)
since the third spatial index is the same as the second.

these signatures are:
    (0)=0, (1)=1  for 2RDM
    (0, 0)=0, (1, 0)=1, (0, 1)=2  for 3RDM
    (0, 0, 0)=0, (1, 0, 0)=1, (0, 1, 0)=2, (0, 0, 1)=4, (1, 0, 1)=5  for 4RDM
'''

def spatsigint_to_unique(i, rank):
    # returns None if invalid
    unique = [0]
    last_t = 0
    for t in int_to_tup(i, rank-1):
        if t and last_t: return None
        if t: unique.append(unique[-1])
        else: unique.append(1+unique[-1])
        last_t = t
    return unique


def generate_spatsigints(rank):
    out = []
    for i in range(2**(rank-1)):
        if spatsigint_to_unique(i, rank) is not None: out.append(i)
    return out

def parity(inds, t_sort=False):
    switches = 0
    last = 0
    inds = list(inds)
    while True:
        for i in range(len(inds)-1):
            if inds[i]>inds[i+1]:
                inds[i], inds[i+1] = inds[i+1], inds[i]
                switches+=1
        if switches==last:
            break
        last = switches
    if t_sort:
        return 1-2*(switches%2), tuple(inds)
    else:
        return 1-2*(switches%2)

def reorder(tup, inds):
    return tuple(tup[i] for i in inds)

def tup_eq(t1, t2):
    return all(i==j for i, j in zip(t1, t2))


def get_reordering_factors(hole_uniques, elec_uniques, hole_spinsig, elec_spinsig):
    '''
    d[unique spatial index ordering] = {
        'hole_reorder': ..., 
        'elec_reorder': ..., 
        'factor': ..., (including parity)
    }
    '''
    d = {}
    '''
    uniques are tuples of length RDM rank which store the similarity
    structure of the spinned element
    e.g. (in NECI FORTRAN indices)
    the spinned element (1, 2, 6, 7) -> (2, 3, 4, 5)
    would have spatials (1, 1, 3, 4)    (1, 2, 2, 3)
    and spinsigs        (1, 0, 0, 1)    (0, 1, 0, 1)
    and uniques         (0, 0, 1, 2)    (0, 1, 1, 2)
    e.g. (in NECI FORTRAN indices: a 2RDM example)
    the spinned element (1, 2) -> (4, 9)
    would have spatials (1, 1)    (2, 5)
    and spinsigs        (1, 0)    (0, 1)
    and uniques         (0, 0)    (0, 1)
    the reorderings are [(0, 1), (1, 0)]
    holes -> holes[(0,1)], elecs -> elecs[(0,1)]:
        spinsigs (1, 0), (0, 1) don't match so no contribution
    holes -> holes[(0,1)], elecs -> elecs[(1,0)]:
        spinsigs (1, 0), (1, 0) do match
        d[(0,0)+(1,0)]['factor'] = -1
    holes -> holes[(1,0)], elecs -> elecs[(0,1)]:
        spinsigs (0, 1), (0, 1) do match
        d[(0,0)+(0,1)]['factor'] = -1
    holes -> holes[(1,0)], elecs -> elecs[(1,0)]:
        spinsigs (0, 1), (1, 0) don't match so no contribution
    '''
    rank = len(hole_uniques)
    assert all(len(x)==rank for x in (elec_uniques, hole_spinsig, elec_spinsig))
    for hole_reorder in itertools.permutations(range(rank)):
        r_hole_uniques = reorder(hole_uniques, hole_reorder)
        r_hole_spinsig = reorder(hole_spinsig, hole_reorder)
        for elec_reorder in itertools.permutations(range(rank)):
            r_elec_uniques = reorder(elec_uniques, elec_reorder)
            r_elec_spinsig = reorder(elec_spinsig, elec_reorder)
            if not tup_eq(r_hole_spinsig, r_elec_spinsig): continue
            p = parity(hole_reorder)*parity(elec_reorder)
            key = r_hole_uniques+r_elec_uniques
            if key not in d: 
                d[key] = {
                    'hole_reorder' : hole_reorder,
                    'elec_reorder' : elec_reorder,
                    'factor' : p
                }
            else: d[key]['factor']+=p
    return d.values()

def is_valid_unique(tup):
    if tup[0]: return False
    if min(tup): return False
    if len(tup)==2: return True
    prev_2, prev_1 = tup[0], tup[1]
    if prev_1-prev_2>1: return False
    for i in range(2, len(tup)):
        if tup[i]==prev_2 and prev_2==prev_1: return False
        prev_2, prev_1 = prev_1, tup[i]
        if prev_1-prev_2>1: return False
    return True

def is_valid_spinned_element(uniques, spinsig):
    rank = len(uniques)
    prev = 2*uniques[0]+spinsig[0]
    for i in range(1, rank):
        this = 2*uniques[i]+spinsig[i]
        if this==prev: return False
    return True

def is_valid(hole_uniques, elec_uniques, hole_spinsig, elec_spinsig):
    # does all validation
    # make sure spin conserving
    if not sum(hole_spinsig)==sum(elec_spinsig): return False
    if not is_valid_unique(hole_uniques): return False
    if not is_valid_unique(elec_uniques): return False
    if not is_valid_spinned_element(hole_uniques, hole_spinsig): return False
    if not is_valid_spinned_element(elec_uniques, elec_spinsig): return False
    return True

for rank in (2, 3):
    cases = []
    spatsigints = generate_spatsigints(rank)
    for h_spatsigint in spatsigints:
        h_uniques = spatsigint_to_unique(h_spatsigint, rank)
        for e_spatsigint in spatsigints:
            e_uniques = spatsigint_to_unique(e_spatsigint, rank)
            for i_spinsigint, (h_spinsigint, e_spinsigint) in enumerate(SpinSigPairs(rank)):
                h_spinsig = int_to_tup(h_spinsigint, rank)
                e_spinsig = int_to_tup(e_spinsigint, rank)
                # reject the invalid combinations
                if not is_valid(h_uniques, e_uniques, h_spinsig, e_spinsig): continue
                # these are the valid combinations, so add them to the relevant arrays
                case = get_reordering_factors(h_uniques, e_uniques, h_spinsig, e_spinsig)
                cases.append((h_spatsigint, e_spatsigint, i_spinsigint, case))

    case_map = -np.ones((max(spatsigints)+1, max(spatsigints)+1, len(SpinSigPairs(rank))), dtype=int)
    perm_end_offsets = []
    factors = []
    hole_perms = []
    elec_perms = []

    for icase, (hspat, espat, ispin, case) in enumerate(cases):
        if not len(perm_end_offsets): perm_end_offsets.append(len(case))
        else: perm_end_offsets.append(len(case)+perm_end_offsets[-1])
        case_map[hspat, espat, ispin] = icase
        for perm in case:
            factors.append(perm['factor'])
            hole_perms.extend(perm['hole_reorder'])
            elec_perms.extend(perm['elec_reorder'])

    nperm = perm_end_offsets[-1]

    perm_end_offsets = np.array(perm_end_offsets)
    factors = np.array(factors)
    hole_perms = np.array(hole_perms).reshape((rank, nperm))
    elec_perms = np.array(elec_perms).reshape((rank, nperm))

    items.append({
        'name': f'c_case_map_{rank}', 
        'data': case_map.transpose((2, 1, 0)), 
        'signed': False, 
        'comment': 'maps the pair-spinsig, the annihilation and creation spatsigs to a spin tracing case id',
        'dimnames': [f'c_npair_spinsig_{rank}', f'c_nspatsig_{rank}', f'c_nspatsig_{rank}']})
    
    items.append({
        'name': f'c_perm_end_offsets_{rank}', 
        'data': perm_end_offsets,
        'signed': False, 
        'comment': f'the end-points on the c_nperm_{rank} dimension for each spin tracing case id',
        'dimnames': [f'c_ncase_{rank}']})

    items.append({
        'name': f'c_factors_{rank}', 
        'data': factors,
        'signed': True, 
        'comment': f'the scalar factors associated with each permutation',
        'dimnames': [f'c_nperm_{rank}']})

    items.append({
        'name': f'c_cre_perms_{rank}', 
        'data': hole_perms.T,
        'signed': False, 
        'comment': f'reordering of the creation indices for the given permutation',
        'dimnames': [f'c_nperm_{rank}', f'c_rank_{rank}']})

    items.append({
        'name': f'c_ann_perms_{rank}', 
        'data': elec_perms.T,
        'signed': False, 
        'comment': f'reordering of the annihilation indices for the given permutation',
        'dimnames': [f'c_nperm_{rank}', f'c_rank_{rank}']})

fname = 'SpinFreeRdmArrays.h'
to_file('spin_free_rdm_arrays', items, [], fname)
assert 0







for rank in (2, 3):
    cases = []
    #for ihole_uniq
    for hole_uniques in itertools.combinations_with_replacement(range(rank), rank):
        int_hole_
        for elec_uniques in itertools.combinations_with_replacement(range(rank), rank):
            for spinsig_icase, (int_hole_spinsig, int_elec_spinsig) in enumerate(SpinSigPairs(rank)):
                hole_spinsig = int_to_tup(int_hole_spinsig, rank)
                elec_spinsig = int_to_tup(int_elec_spinsig, rank)
                # reject the invalid combinations
                if not is_valid(hole_uniques, elec_uniques, hole_spinsig, elec_spinsig): continue
                cases.append(get_reordering_factors(hole_uniques, elec_uniques, hole_spinsig, elec_spinsig))

    perm_end_offsets = ParameterArray()
    for case in cases:
        if not len(perm_end_offsets):
            perm_end_offsets.append(len(case))
        else:
            perm_end_offsets.append(len(case)+perm_end_offsets[-1])

    nperm = perm_end_offsets[-1]

    hole_perms = ParameterArray()
    elec_perms = ParameterArray()
    factors = ParameterArray()
    for case in cases:
        for perm in case:
            factors.append(perm['factor'])
            hole_perms.extend(perm['hole_reorder'])
            elec_perms.extend(perm['elec_reorder'])
    hole_perms.set_shape(rank, nperm)
    elec_perms.set_shape(rank, nperm)

    pm.add(perm_end_offsets, 'perm_end_offsets_{}rdm'.format(rank))
    pm.add(factors, 'factors_{}rdm'.format(rank))
    pm.add(hole_perms, 'hole_perms_{}rdm'.format(rank))
    pm.add(elec_perms, 'elec_perms_{}rdm'.format(rank))

pm.add_comment(\
'''{}
this file was generated by meta/{} and should not be modified by hand.
'''.format('robert.anderson@kcl.ac.uk', __file__))
pm.render('..')
