import itertools

'''
ket and bra are lists of spin orbital indices in ascending order
the ket is turned into the bra by the premultiplication of an excitation operator
this operator has nremoved annihilation operators, and ninserted creation operators
due to fermion antisymmetry, this process will return the bra multiplied by an integer
power (nperm) of -1. where nperm is the number of antisymmetric exchanges required to
eliminate the inner product.
'''
def half_parity(det, removed):
    nperm = 0
    nocc_passed = 0
    for occ in det:
        if occ in removed:
            nperm+=nocc_passed
            if occ==removed[-1]: break
        else:
            nocc_passed+=1
    nperm+=sum(i>max(det) for i in removed)*nocc_passed
    return (-1.0)**nperm

def parity(ket, bra):
    assert ket==sorted(ket)
    assert bra==sorted(bra)
    removed = sorted(set(ket)-set(bra))
    inserted = sorted(set(bra)-set(ket))
    return half_parity(ket, removed)*half_parity(bra, inserted)

norb = 8
with open('parity_{}.txt'.format(norb), 'w') as f:
    for ket in itertools.product(*[[0,1]]*norb):
        lket = list(filter(lambda i:ket[i], range(norb)))
        for bra in itertools.product(*[[0,1]]*norb):
            lbra = list(filter(lambda i:bra[i], range(norb)))
            assert parity(lket, lbra)==parity(lbra, lket)
            f.write('{} {} {}\n'.format(
                parity(lket, lbra),
                ' '.join(map(str, ket)),
                ' '.join(map(str, bra))))
