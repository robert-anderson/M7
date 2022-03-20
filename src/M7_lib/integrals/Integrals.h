//
// Created by Robert John Anderson on 2020-02-01.
//

#ifndef M7_INTEGRALS_H
#define M7_INTEGRALS_H

static constexpr size_t trig(size_t i, size_t j){
    return i+(j*(j+1))/2;
}
class Integrals {
public:
    const size_t m_nsite; // the number of spatial orbitals
    const bool m_spin_res; // are we storing integrals with spin-resolved indices?
    const size_t m_nspinorb; // the number of spin orbitals
    const size_t m_nintind; // the extent of the stored integral indices
protected:
    Integrals(size_t nsite, bool spin_res) :
            m_nsite(nsite), m_spin_res(spin_res),
            m_nspinorb(2*nsite), m_nintind(spin_res ? m_nspinorb : m_nsite){}

    size_t spinorb(const size_t &ispat, const size_t &ispin) const {
        if (m_spin_res) return ispin * m_nsite + ispat;
        return ispat;
    }
};

#endif //M7_INTEGRALS_H
