import numpy as np
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

h5_file = h5py.File('M7.h5', 'r')

h5_rdms = h5_file['archive']['rdms']
with open('exact_rdms.pkl', 'rb') as f: py_rdms = pkl.load(f)

print(h5_rdms.keys())
print(py_rdms.keys())

for key in ('sf_3300',):
    h5_rdm = load_hdf5_rdm(h5_rdms[key])
    py_rdm = py_rdms[key]
    print(np.einsum('ijkijk->', h5_rdm))
    print(np.einsum('ijkijk->', py_rdm))
    max_diff = max(np.abs(h5_rdm-py_rdm).flatten())
    print(key, max_diff)
