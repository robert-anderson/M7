from pyscf import M, mrpt, mcscf, scf, fci, symm, tools
import sys, os, h5py
sys.path.append(os.path.abspath('.'))
import numpy as np
from functools import reduce

d = {
    'mol_kwargs': {
        'atom':
        '''
        H 0.000000 0.000000 0.000000
        F 0.000000 0.000000 1.4
        ''',
        'basis': {'H': 'cc-pvdz', 'F': 'cc-pvdz'},
        'symmetry': True,
        'symmetry_subgroup': 'C2v',
        'verbose': 4,
        'charge': 0,
        'spin': 0
        }
    }

def fcidump(mc, mol, filename='FCIDUMP.new'):
    IRREP_MAP = {
        'D2h': (1, 4, 6, 7, 8, 5, 3, 2),
        'C2v': (1, 4, 2, 3),
        'C2h': (1, 4, 2, 3),
        'D2' : (1, 4, 3, 2),
        'Cs' : (1, 2),
        'C2' : (1, 2),
        'Ci' : (1, 2),
        'C1' : (1,)
    }
    eri_cas = mc.get_h2eff(mc.mo_coeff)
    h1eff, energy_core = mc.get_h1eff(mc.mo_coeff)
    orbsym = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, mc.mo_coeff)
    orbsym = [symm.param.IRREP_ID_TABLE[mol.groupname][i] for i in orbsym]
    orbsym = [IRREP_MAP[mol.groupname][i] for i in orbsym]
    orbsym = orbsym[mc.ncore:mc.ncore+h1eff.shape[0]]
    tools.fcidump.from_integrals(filename, h1eff, eri_cas, h1eff.shape[0], mc.nelecas, energy_core, mol.ms, orbsym)

'''
output Fock matrix in Molcas HDF5 format
'''
def fockdump(fock, fname):
    inds = []
    values = []
    for i in range(fock.shape[0]):
        for j in range(0, i+1):
            if fock[i, j]:
                inds.append((i+1, j+1))
                values.append(fock[i,j])
    inds = np.array(inds)
    values = np.array(values)

    f = h5py.File(fname, 'w')
    f['ACT_FOCK_INDEX'] = inds
    f['ACT_FOCK_VALUES'] = values

nelecas = 6
ncas = 5

mol = M(**d['mol_kwargs'])
nelecb = (nelecas-mol.spin)//2
neleca = nelecas - nelecb

myhf = scf.ROHF(mol)
myhf.kernel()
myhf.analyze()

mc = mcscf.CASCI(myhf, ncas, nelecas)
mc.kernel()

ao_rdm = mc.make_rdm1()[mc.ncore:mc.ncore+ncas, mc.ncore:mc.ncore+ncas]
mo_rdm = mc.fcisolver.make_rdm1(mc.fcisolver.ci, ncas, nelecas)

fock_ao = mc.get_fock(mc.mo_coeff, mc.fcisolver.ci, mc.get_h2eff(mc.mo_coeff), mo_rdm)
fock_mo = reduce(np.dot, (mc.mo_coeff.T, fock_ao, mc.mo_coeff))
cas_fock_mo = fock_mo[mc.ncore:mc.ncore+ncas, mc.ncore:mc.ncore+ncas]

dm1, dm2, dm3, dm4 = fci.rdm.make_dm1234('FCI4pdm_kern_sf', mc.fcisolver.ci, mc.fcisolver.ci, ncas, (neleca, nelecb))
dm1, dm2, dm3, dm4 = fci.rdm.reorder_dm1234(dm1, dm2, dm3, dm4)
dm4f = np.einsum('iajbkcld,ld->iajbkc', dm4, cas_fock_mo)

assert ao_rdm.shape==mo_rdm.shape

fcidump(mc, mol, 'FCIDUMP')
fockdump(cas_fock_mo, 'fock.h5')

# transpose higher-rank arrays such that the creation inds are first, then the annihilations
dms = {
    'sf_1100' : np.array(dm1),
    'sf_2200' : np.array(dm2).transpose(0, 2, 1, 3),
    'sf_3300' : np.array(dm3).transpose(0, 2, 4, 1, 3, 5),
    'sf_4400f': np.array(dm4f).transpose(0, 2, 4, 1, 3, 5),
}

import pickle as pkl
with open('exact_rdms.pkl', 'wb') as f:
    pkl.dump(dms, f)
