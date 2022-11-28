from pyscf import gto, scf, fci, ao2mo
import numpy as np
from functools import reduce

#
# 1. open-shell system
#
mol = gto.M(
    atom = '''N        0.000000    0.000000    0.000000''',
    basis = '6-31g',
    charge = 0,
    spin = 3  # = 2S = spin_up - spin_down
)

mf = scf.UHF(mol)
mf.kernel()

print(mf.mo_coeff.shape)
cisolver = fci.FCI(mf)
cisolver.nroots = 5
efci = cisolver.kernel()[0]
for i in range(cisolver.nroots):
    print('E(UHF-FCI) = %.12f' % efci[i])


''' 
first prepare the spin-resolved, spin-minor FCIDUMP
'''

fname = 'FCIDUMP.spin_minor'

orbs = mf.mo_coeff
nmo = orbs.shape[-1]
orbsym = [1,]*nmo
nelec = mol.nelectron
tol = 1e-12
ms = mol.spin
eri_aaaa = ao2mo.restore(8,ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[0],orbs[0]), compact=False),nmo)
eri_bbbb = ao2mo.restore(8,ao2mo.incore.general(mf._eri, (orbs[1],orbs[1],orbs[1],orbs[1]), compact=False),nmo)
eri_aabb = ao2mo.restore(8,ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[1],orbs[1]), compact=False),nmo)
eri_bbaa = ao2mo.restore(8,ao2mo.incore.general(mf._eri, (orbs[1],orbs[1],orbs[0],orbs[0]), compact=False),nmo)
h_core = mf.get_hcore(mol)
h_aa = reduce(np.dot, (orbs[0].T, h_core, orbs[0]))
h_bb = reduce(np.dot, (orbs[1].T, h_core, orbs[1]))
nuc = mol.energy_nuc()
float_format = ' %.16g'

# Stupidly, NECI wants its orbitals as a,b,a,b,a,b rather than aaaabbbb
# Reorder things so this is the case
#assert(len(orbsym) % 2 == 0)
#orbsym_reorder = [i for tup in zip(orbsym[:len(orbsym)//2], orbsym[len(orbsym)//2:]) for i in tup]
orbsym_reorder = orbsym
a_inds = [i*2+1 for i in range(orbs[0].shape[1])]
b_inds = [i*2+2 for i in range(orbs[1].shape[1])]


eri_inds = []
ij = 0
for i in range(nmo):
    for j in range(0, i+1):
        kl = 0
        for k in range(0, i+1):
            for l in range(0, k+1):
                if ij >= kl:
                    eri_inds.append((i, j, k, l))
                kl += 1
        ij += 1

h_inds = []
for i in range(nmo):
    for j in range(0, i+1):
        h_inds.append((i, j))


with open(fname, 'w') as fout:
    fout.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' % (nmo, nelec, ms))
    if orbsym is not None and len(orbsym_reorder) > 0:
        fout.write('  ORBSYM=%s\n' % ','.join([str(x) for x in orbsym_reorder]))
    else:
        fout.write('  ORBSYM=%s\n' % ('1,' * 2*nmo))
    fout.write('  ISYM=1, UHF=TRUE\n')
    fout.write(' &END\n')
    # Assume 8-fold symmetry
    output_format = float_format + ' %4d %4d %4d %4d\n'
    for ijkl, (i, j, k, l) in enumerate(eri_inds):
        if abs(eri_aaaa[ijkl]) > tol:
            fout.write(output_format % (eri_aaaa[ijkl], a_inds[i], a_inds[j], a_inds[k], a_inds[l]))
        if abs(eri_bbbb[ijkl]) > tol:
            fout.write(output_format % (eri_bbbb[ijkl], b_inds[i], b_inds[j], b_inds[k], b_inds[l]))
        if abs(eri_aabb[ijkl]) > tol:
            fout.write(output_format % (eri_aabb[ijkl], a_inds[i], a_inds[j], b_inds[k], b_inds[l]))
        if abs(eri_bbaa[ijkl]) > tol:
            fout.write(output_format % (eri_bbaa[ijkl], b_inds[i], b_inds[j], a_inds[k], a_inds[l]))
    h_aa = h_aa.reshape(nmo,nmo)
    h_bb = h_bb.reshape(nmo,nmo)
    output_format = float_format + ' %4d %4d  0  0\n'
    for ij, (i, j) in enumerate(h_inds):
        if abs(h_aa[i,j]) > tol:
            fout.write(output_format % (h_aa[i,j], a_inds[i], a_inds[j]))
        if abs(h_bb[i,j]) > tol:
            fout.write(output_format % (h_bb[i,j], b_inds[i], b_inds[j]))
    output_format = float_format + '  0  0  0  0\n'
    fout.write(output_format % nuc)

fname = 'FCIDUMP.spin_major'
'''
then spin major
'''
a_inds = [i + 1 for i in range(orbs[0].shape[1])]
b_inds = [nmo + i + 1 for i in range(orbs[1].shape[1])]
with open(fname, 'w') as fout:
    fout.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' % (nmo, nelec, ms))
    if orbsym is not None and len(orbsym_reorder) > 0:
        fout.write('  ORBSYM=%s\n' % ','.join([str(x) for x in orbsym_reorder]))
    else:
        fout.write('  ORBSYM=%s\n' % ('1,' * 2*nmo))
    fout.write('  ISYM=1, UHF=TRUE\n')
    fout.write(' &END\n')
    # Assume 8-fold symmetry
    output_format = float_format + ' %4d %4d %4d %4d\n'
    for ijkl, (i, j, k, l) in enumerate(eri_inds):
        if abs(eri_aaaa[ijkl]) > tol:
            fout.write(output_format % (eri_aaaa[ijkl], a_inds[i], a_inds[j], a_inds[k], a_inds[l]))
        if abs(eri_bbbb[ijkl]) > tol:
            fout.write(output_format % (eri_bbbb[ijkl], b_inds[i], b_inds[j], b_inds[k], b_inds[l]))
        if abs(eri_aabb[ijkl]) > tol:
            fout.write(output_format % (eri_aabb[ijkl], a_inds[i], a_inds[j], b_inds[k], b_inds[l]))
        if abs(eri_bbaa[ijkl]) > tol:
            fout.write(output_format % (eri_bbaa[ijkl], b_inds[i], b_inds[j], a_inds[k], a_inds[l]))
    h_aa = h_aa.reshape(nmo,nmo)
    h_bb = h_bb.reshape(nmo,nmo)
    output_format = float_format + ' %4d %4d  0  0\n'
    for ij, (i, j) in enumerate(h_inds):
        if abs(h_aa[i,j]) > tol:
            fout.write(output_format % (h_aa[i,j], a_inds[i], a_inds[j]))
        if abs(h_bb[i,j]) > tol:
            fout.write(output_format % (h_bb[i,j], b_inds[i], b_inds[j]))
    output_format = float_format + '  0  0  0  0\n'
    fout.write(output_format % nuc)


''' 
next prepare the Molpro block formatted FCIDUMP
'''
fname = 'FCIDUMP.spin_blocks'

output_format = float_format + '  0  0  0  0\n'
delimit_line = output_format % 0.0

with open(fname, 'w') as fout:
    fout.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' % (nmo, nelec, ms))
    if orbsym is not None and len(orbsym_reorder) > 0:
        fout.write('  ORBSYM=%s\n' % ','.join([str(x) for x in orbsym_reorder]))
    else:
        fout.write('  ORBSYM=%s\n' % ('1,' * 2*nmo))
    fout.write('  ISYM=1, UHF=TRUE\n')
    fout.write(' &END\n')
    # Assume 8-fold symmetry
    output_format = float_format + ' %4d %4d %4d %4d\n'
    for ijkl, (i, j, k, l) in enumerate(eri_inds):
        if abs(eri_aaaa[ijkl]) > tol:
            fout.write(output_format % (eri_aaaa[ijkl], i+1, j+1, k+1, l+1))
    fout.write(delimit_line)

    for ijkl, (i, j, k, l) in enumerate(eri_inds):
        if abs(eri_aabb[ijkl]) > tol:
            fout.write(output_format % (eri_aabb[ijkl], i+1, j+1, k+1, l+1))
        if (i,j)!=(k,l):
            if abs(eri_bbaa[ijkl]) > tol:
                fout.write(output_format % (eri_bbaa[ijkl], k+1, l+1, i+1, j+1))
    fout.write(delimit_line)

    for ijkl, (i, j, k, l) in enumerate(eri_inds):
        if abs(eri_bbbb[ijkl]) > tol:
            fout.write(output_format % (eri_bbbb[ijkl], i+1, j+1, k+1, l+1))
    fout.write(delimit_line)

    h_aa = h_aa.reshape(nmo,nmo)
    h_bb = h_bb.reshape(nmo,nmo)
    output_format = float_format + ' %4d %4d  0  0\n'
    for ij, (i, j) in enumerate(h_inds):
        if abs(h_aa[i,j]) > tol:
            fout.write(output_format % (h_aa[i,j], i+1, j+1))
    fout.write(delimit_line)
    for ij, (i, j) in enumerate(h_inds):
        if abs(h_bb[i,j]) > tol:
            fout.write(output_format % (h_bb[i,j], i+1, j+1))
    fout.write(delimit_line)
    output_format = float_format + '  0  0  0  0\n'
    fout.write(output_format % nuc)




