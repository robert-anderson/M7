//
// Created by Robert John Anderson on 2020-02-01.
//

#ifndef M7_INTEGRALS_H
#define M7_INTEGRALS_H

constexpr size_t trig(size_t i, size_t j){
    return i+(j*(j+1))/2;
}

class Integrals {
public:
    const size_t m_norb; // the extent of the stored integral indices
    const bool m_spin_resolved; // are we storing integrals with spin-resolved indices?
    const size_t m_nspinorb; // the number of spin orbitals
    const size_t m_nspatorb; // the number of spatial orbitals
protected:
    Integrals(size_t norb, bool spin_resolved) :
    m_norb(norb), m_spin_resolved(spin_resolved),
    m_nspinorb(spin_resolved?norb:2*norb),
    m_nspatorb(spin_resolved?norb/2:norb){}

    size_t spinorb(const size_t &ispat, const size_t &ispin) const {
        if (m_spin_resolved) return ispin*m_nspatorb+ispat;
        return ispat;
    }
};


#endif //M7_INTEGRALS_H
