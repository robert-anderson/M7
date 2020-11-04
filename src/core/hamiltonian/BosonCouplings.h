//
// Created by RJA on 22/09/2020.
//

#ifndef M7_BOSONCOUPLINGS_H
#define M7_BOSONCOUPLINGS_H

#if 0

#include <src/core/basis/PermanentConnection.h>
#include "Hamiltonian.h"


struct FermionDensityCoupling {
    virtual defs::ham_t operator()(const size_t &p, const size_t &q, const size_t &n) const = 0;
};

struct FermionDensityCouplingArray {
    // TODO
};

struct FermionDensityCouplingParameter : FermionDensityCoupling{
    defs::ham_t m_v;
    FermionDensityCouplingParameter(defs::ham_t v) : m_v(v) {}
    defs::ham_t operator()(const size_t &p, const size_t &q, const size_t &n) const override {
        return m_v;
    }
};


struct BosonModeFrequencies {
    virtual defs::ham_t operator()(const size_t &n) const = 0;
};

struct BosonModeFrequenciesArray {
    // TODO
};

struct BosonModeFrequenciesParameter : BosonModeFrequencies{
    defs::ham_t m_omega;
    BosonModeFrequenciesParameter(defs::ham_t omega) : m_omega(omega) {}
    defs::ham_t operator()(const size_t &n) const override {
        return m_omega;
    }
};


class BosonCouplings {

    const size_t m_nocc_cutoff;
    const size_t m_nmode;
    Hamiltonian *m_ham;
    std::unique_ptr<FermionDensityCoupling> m_v;
    std::unique_ptr<BosonModeFrequencies> m_omega;

public:
    BosonCouplings(size_t nocc_cutoff, size_t nmode, defs::ham_t v, defs::ham_t omega) :
            m_nocc_cutoff(nocc_cutoff), m_nmode(nmode),
            m_v(std::unique_ptr<FermionDensityCoupling>(new FermionDensityCouplingParameter(v))),
            m_omega(std::unique_ptr<BosonModeFrequencies>(new BosonModeFrequenciesParameter(omega))){}

    defs::ham_t v(const size_t& p, const size_t& q, const size_t& imode) const {
        return (*m_v)(p, q, imode);
    }

    defs::ham_t omega(const size_t& imode) const {
        return (*m_omega)(imode);
    }

    const size_t& nocc_cutoff() const {
        return m_nocc_cutoff;
    }
    const size_t& nmode() const {
        return m_nmode;
    }
    defs::ham_t get_element_0(const PermanentConnection &permconn) const {
        defs::ham_t res = 0;
        for (size_t imode=0ul; imode<m_nmode; ++imode)
            res+= omega(imode)*permconn.com(imode);
        return res;
    }

    defs::ham_t get_element_1(const AntisymConnection &detconn, const PermanentConnection &permconn) const {
        const auto n = permconn.changed_modes(0);
        // bosons don't couple to higher fermion excitations (yet?)
        switch (detconn.nexcit()) {
            case 0: {
                defs::ham_t res = 0;
                for (size_t iocc = 0ul; iocc < detconn.ncom(); ++iocc) {
                    auto p = detconn.com()[iocc];
                    res += v(p, p, n);
                }
                return res;
            }
            case 1:{
                auto p = detconn.cre()[0];
                auto q = detconn.ann()[0];
                return v(p, q, n);
            }
            default:
                return 0;
        }
    }

    defs::ham_t get_element(const AntisymConnection &detconn,
                            const PermanentConnection &permconn) const {

        switch (permconn.nchanged_mode()) {
            case 0:
                return get_element_0(permconn);
            case 1:
                return get_element_1(detconn, permconn);
            default:
                return 0;
        }
    }

};


#endif //M7_BOSONCOUPLINGS_H
#endif //M7_BOSONCOUPLINGS_H
