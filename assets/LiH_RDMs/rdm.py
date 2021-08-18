from pyscf import M, mrpt, mcscf, scf, fci, symm, fciqmcscf, tools, ao2mo
import sys, os
from functools import reduce
sys.path.append(os.path.abspath('.'))
import numpy

d = {
        'mol_kwargs': {
            'atom':
            '''
            H 0.000000 0.000000 0.000000
            Li 0.000000 0.000000 1.4
            ''',
            'basis': {'H': 'sto-3g', 'Li': 'sto-3g'},
            'symmetry': True,
            'symmetry_subgroup': 'C2v',
            'verbose': 4,
            'charge': 0,
            'spin': 0
        }
    }

mol = M(**d['mol_kwargs'])
myhf = scf.RHF(mol)
myhf.kernel()
myhf.analyze()

myfci = fci.FCI(myhf)
fcie, fcivec = myfci.kernel()

nelecb = (mol.nelectron-mol.spin)//2
neleca = mol.nelectron - nelecb
norb = len(myhf.mo_energy)

dm1, dm2, dm3 = fci.rdm.make_dm123('FCI3pdm_kern_sf', fcivec, fcivec, norb, (neleca, nelecb))

numpy.set_printoptions(precision=4)

print('1RDM diagonal: {}'.format(dm1.diagonal()))
print(dm1)


h1e = reduce(numpy.dot, (myhf.mo_coeff.T, myhf.get_hcore(), myhf.mo_coeff))
eri = mol.intor('int2e').reshape(norb, norb, norb, norb)
eri = ao2mo.kernel(ao2mo.restore(8, eri, norb), myhf.mo_coeff, verbose=0, compact=False)


tools.fcidump.from_integrals('FCIDUMP', h1e, eri, norb, mol.nelectron, myhf.energy_nuc(), 0, [1,]*norb)

dms = {
        '1': numpy.array(dm1),
        '2': numpy.array(dm2),
        '3': numpy.array(dm3)
        }

import pickle as pkl
with open('exact_rdms.pkl', 'wb') as f:
    pkl.dump(dms, f)
