import numpy as np
import sys, itertools
sys.path.append('/home/rja/CLionProjects/M7/python')
from pyscf import fci
from rdms import *

def elementwise_compare(a, b, tol=1e-6):
    rank = len(a.shape)//2
    norb = a.shape[0]
    assert len(a.shape)==len(b.shape)
    for inds in itertools.product(*tuple(range(norb) for i in range(rank*2))):
        if abs(a[inds]-b[inds])>tol:
            print('{}   {:.5f}   {:.5f}     {:.3f}'.format(inds, a[inds], b[inds], a[inds]/b[inds]))


def cosine_diff(a1, a2):
    a1 = a1.flatten()
    a2 = a2.flatten()
    return 1-np.dot(a1, a2)/np.sqrt(np.dot(a1, a1)*np.dot(a2,a2))


import pickle as pkl
exact_fname = '/home/rja/CLionProjects/M7/assets/LiH_RDMs/exact_rdms.pkl'
with open(exact_fname, 'rb') as f:
    exact_rdms = pkl.load(f)


rdm1 = load_spin_resolved_rdm('M7.save.h5', 1)
nelec = sum(rdm1.diagonal())
rdm1_restored = restore_perm_syms(rdm1, True, False)
rdm1_sf = spin_resolved_to_spinfree(rdm1_restored, True)

np.set_printoptions(precision=4)
print("from M7:")
print(rdm1_sf)
print("exact: ")
print(exact_rdms['1'])
assert cosine_diff(rdm1_sf, exact_rdms['1']) < 1e-13

rdm2 = load_spin_resolved_rdm('M7.save.h5', 2)
rdm2_restored = restore_perm_syms(rdm2, True, False)
rdm2_sf = spin_resolved_to_spinfree(rdm2_restored, True)
_, exact_rdm2_sf = reorder_rdm12(exact_rdms['1'], exact_rdms['2'], False)

elementwise_compare(rdm2_sf, exact_rdm2_sf)
assert cosine_diff(rdm2_sf, exact_rdm2_sf) < 1e-13
rdm1_2 = one_from_two_rdm(rdm2_sf, nelec)
assert np.allclose(rdm1_2, rdm1_sf)
