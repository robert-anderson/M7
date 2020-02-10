import numpy as np
import itertools

def file_to_dict(fname, test_eri4fold=None):
    d = dict()
    norb = 0
    with open(fname, 'r') as f:
        for line in f.readlines():
            try: inds = tuple(map(lambda i : int(i)-1, line.split()[-4:]))
            except ValueError: continue
            if any(i<0 for i in inds): continue
            if max(inds)>norb: norb=max(inds)
            # transform to physics notation
            #inds = tuple(inds[i] for i in (0,2,1,3))
            value = line.split()[0][1:-1].split(',')
            value = complex(float(value[0]), float(value[1]))
            d[inds] = value
            if test_eri4fold is not None:
                print inds
                assert abs(test_eri4fold.symm_and_get(inds)-value) < 1e-14

    return d, norb+1

def trig(i, j): return i + (j*(j+1))/2

dump_dict, norb = file_to_dict('DHF_Be_STO-3G/FCIDUMP')

class Eri4Fold:
    def __init__(self, dump_dict, norb):
        self.nslots_8fold = trig(0, trig(0, norb))
        self.nslots_4fold = 2*self.nslots_8fold
        self.data = np.zeros(self.nslots_4fold, dtype=np.complex)
        for inds, value in dump_dict.iteritems():
            self.symm_and_set(inds, value)
        
    def set_value(self, flat_ind, value, tconj, eps=1e-14):
        value = np.conj(value) if tconj else value
        if abs(self.data[flat_ind])> eps:
            assert abs(self.data[flat_ind]-value) < eps, '|{} - {}| = {}'.format(self.data[flat_ind], value, abs(self.data[flat_ind]-value))
        self.data[flat_ind] = value

    def get_value(self, flat_ind, tconj):
        return np.conj(self.data[flat_ind]) if tconj else self.data[flat_ind]

    def symm_and_set(self, inds, value):
        '''
        Symmetries:
        (I) : (ij|kl) = (kl|ij) 
              indistinguishability, invariance to an exchange of dummy indices in the integrand
        (R) : (ij|kl) = (ij|lk)
              real orbitals
        (H) : (ij|kl) = (ji|lk)
              hermiticity

        I and H commute:        IH(ij|kl) = I(ji|lk) = (lk|ji)
                                HI(ij|kl) = H(kl|ij) = (lk|ji)

        I and R do not commute: IR(ij|kl) = I(ij|lk) = (lk|ij)
                                RI(ij|kl) = R(kl|ij) = (kl|ji) = H(lk|ij)

        H and R commute:        HR(ij|kl) = H(ij|lk) = (ji|kl)
                                RH(ij|kl) = R(ji|lk) = (ji|kl)
        '''
        i, j, k, l = inds
        if i<=j:
            if k<=l:
                if trig(i, j) <= trig(k, l):
                    # (ij|kl)
                    self.set_value(trig(trig(i, j), trig(k, l)), value, 0)
                else:
                    # (kl|ij) : (I)
                    self.set_value(trig(trig(k, l), trig(i, j)), value, 0)
            else:
                if trig(i, j) <= trig(l, k):
                    # (ij|lk) : (R)
                    self.set_value(self.nslots_8fold + trig(trig(i, j), trig(l, k)), value, 0)
                else:
                    # (lk|ij) : (IR)
                    self.set_value(self.nslots_8fold + trig(trig(l, k), trig(i, j)), value, 1)
        else:
            if k<=l:
                if trig(j, i) <= trig(k, l):
                    # (ji|kl) : (HR)
                    self.set_value(self.nslots_8fold + trig(trig(j, i), trig(k, l)), value, 1)
                else:
                    # (kl|ji) : (HIR) = (RI)
                    self.set_value(self.nslots_8fold + trig(trig(k, l), trig(j, i)), value, 0)
            else:
                if trig(j, i) <= trig(l, k):
                    # (ji|lk) : (H)
                    self.set_value(trig(trig(j, i), trig(l, k)), value, 1)
                else:
                    # (lk|ji) : (HI)
                    self.set_value(trig(trig(l, k), trig(j, i)), value, 1)

    def symm_and_get(self, inds):
        i, j, k, l = inds
        if i<=j:
            if k<=l:
                if trig(i, j) <= trig(k, l):
                    # (ij|kl)
                    return self.get_value(trig(trig(i, j), trig(k, l)), 0)
                else:
                    # (kl|ij) : (I)
                    return self.get_value(trig(trig(k, l), trig(i, j)), 0)
            else:
                if trig(i, j) <= trig(l, k):
                    # (ij|lk) : (R)
                    return self.get_value(self.nslots_8fold + trig(trig(i, j), trig(l, k)), 0)
                else:
                    # (lk|ij) : (IR)
                    return self.get_value(self.nslots_8fold + trig(trig(l, k), trig(i, j)), 1)
        else:
            if k<=l:
                if trig(j, i) <= trig(k, l):
                    # (ji|kl) : (HR)
                    return self.get_value(self.nslots_8fold + trig(trig(j, i), trig(k, l)), 1)
                else:
                    # (kl|ji) : (HIR) = (RI)
                    return self.get_value(self.nslots_8fold + trig(trig(k, l), trig(j, i)), 0)
            else:
                if trig(j, i) <= trig(l, k):
                    # (ji|lk) : (H)
                    return self.get_value(trig(trig(j, i), trig(l, k)), 1)
                else:
                    # (lk|ji) : (HI)
                    return self.get_value(trig(trig(l, k), trig(j, i)), 1)
        assert 0

eri4fold = Eri4Fold(dump_dict, norb)
file_to_dict('DHF_Be_STO-3G/FCIDUMP', eri4fold)


