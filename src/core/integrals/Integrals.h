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


struct IntegralArray {
    static constexpr size_t trig(size_t i, size_t j){
        return i+(j*(j+1))/2;
    }
    static constexpr size_t nintind(size_t nsite, bool spin_res){
        return spin_res ? 2*nsite : nsite;
    }
    /**
     * number of spatial orbitals in the system
     */
    const size_t m_nsite;
    /**
     * "spin resolved" true if Sz symmetry is not conserved in the one-particle basis e.g. UHF, DiracHF
     */
    const bool m_spin_res;
    /**
     * number of spin orbitals
     */
    const size_t m_nspinorb;
    /**
     * extent of the stored integral indices
     */
    const size_t m_nintind;
protected:
    IntegralArray(size_t nsite, bool spin_res) :
        m_nsite(nsite), m_spin_res(spin_res),
            m_nspinorb(2*nsite), m_nintind(spin_res ? m_nspinorb : m_nsite){}

    size_t spinorb(const size_t &ispat, const size_t &ispin) const {
        if (m_spin_res) return ispin * m_nsite + ispat;
        return ispat;
    }
};

struct IntegralArray1e : IntegralArray {
    std::vector<defs::ham_t> m_data;
    IntegralArray1e(size_t nsite, bool spin_res, size_t size):
        IntegralArray(nsite, spin_res), m_data(size, {}){}
    virtual void get(const size_t& i, const size_t& j, defs::ham_comp_t& elem) const {}
    virtual void get(const size_t& i, const size_t& j, std::complex<defs::ham_comp_t>& elem) const {}
    defs::ham_t get(const size_t& i, const size_t& j) const {
        defs::ham_t elem;
        get(i, j, elem);
        return elem;
    }

    defs::ham_t& operator[](const size_t& iflat) {
        DEBUG_ASSERT_LT(iflat, m_data.size(), "integral index OOB");
        return m_data[iflat];
    }

    const defs::ham_t& operator[](const size_t& iflat) const {
        DEBUG_ASSERT_LT(iflat, m_data.size(), "integral index OOB");
        return m_data[iflat];
    }
};

struct IntegralArray1e_1fold : IntegralArray1e {
    IntegralArray1e_1fold(size_t nsite, bool spin_res):
        IntegralArray1e(nsite, spin_res, utils::pow<2>(nintind(nsite, spin_res))){}
    void get(const size_t &i, const size_t &j, defs::ham_comp_t &elem) const override {
        elem = (*this)[i*m_nsite+j];
    }
    void get(const size_t &i, const size_t &j, std::complex<defs::ham_comp_t> &elem) const override {
        elem = (*this)[i*m_nsite+j];
    }
};

struct IntegralArray1e_2fold : IntegralArray1e {
    IntegralArray1e_2fold(size_t nsite, bool spin_res):
        IntegralArray1e(nsite, spin_res, IntegralArray::trig(0, nintind(nsite, spin_res))){}

    void get(const size_t &i, const size_t &j, defs::ham_comp_t &elem) const override {
        elem = i<=j ? (*this)[trig(i, j)] : (*this)[trig(j, i)];
    }
    void get(const size_t &i, const size_t &j, std::complex<defs::ham_comp_t> &elem) const override {
        elem = i<=j ? (*this)[trig(i, j)] : std::conj((*this)[trig(j, i)]);
    }
};

//struct IntegralArray2e : IntegralArray {};
//
//struct IntegralArray2e_8fold : IntegralArray2e {};


#endif //M7_INTEGRALS_H
