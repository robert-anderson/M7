from pyscf import M, mrpt, mcscf, scf, fci, symm, fciqmcscf
import sys, os
sys.path.append(os.path.abspath('.'))
import numpy as np

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


dm1, dm2, dm3 = fci.rdm.make_dm123('FCI3pdm_kern_sf', mc.fcisolver.ci, mc.fcisolver.ci, ncas, (neleca, nelecb))
dm1, dm2, dm3 = fci.rdm.reorder_dm123(dm1, dm2, dm3)

assert ao_rdm.shape==mo_rdm.shape
print(mo_rdm.diagonal())
print(sum(mo_rdm.diagonal()))
print(mo_rdm)

dms = {
        '1': np.array(dm1),
        '2': np.array(dm2),
        '3': np.array(dm3)
        }

import pickle as pkl
with open('exact_rdms.pkl', 'wb') as f:
    pkl.dump(dms, f)

'''
mc = mcscf.CASCI(myhf, ncas, nelecas)
mc.fcisolver = fciqmcscf.fciqmc.FCIQMCCI(mol)
mc.kernel()
'''
