from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run(assets=['HF_RDMs/FCIDUMP', 'HF_RDMs/fock_diag.h5', 'HF_RDMs/exact_rdms.pkl'])
check_nw()
check_rdm_archives()

import h5py
import pickle as pkl
def load_hdf5_rdm(group):
    inds = np.array(group['indices'])
    values = np.array(group['values'])
    extent = max(inds.flatten())+1
    nind = inds.shape[1]
    rdm = np.zeros((extent,)*nind)
    for i, row in enumerate(inds): rdm[tuple(row)] = values[i]
    return rdm

h5_fname = make_local_name('M7.rdm.h5')
pkl_fname = make_local_name('exact_rdms.pkl')
h5_file = h5py.File(h5_fname, 'r')

h5_rdms = h5_file['spinfree']
with open(pkl_fname, 'rb') as f: py_rdms = pkl.load(f)

h5_rdm = load_hdf5_rdm(h5_rdms['4400f'])
py_rdm = py_rdms['4400fd']
diff = np.abs(h5_rdm-py_rdm).flatten()
imax = np.argmax(diff)
if diff[imax] > 1e-4: fail(f'F*4RDM does not match exact value {h5_rdm.flatten()[imax], py_rdm.flatten()[imax]}')
