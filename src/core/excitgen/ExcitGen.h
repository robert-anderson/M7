//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITGEN_H
#define M7_EXCITGEN_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/connection/Connections.h>
#include "src/core/sample/PRNG.h"
#include "CachedOrbs.h"

struct ExcitGen {

protected:
    const Hamiltonian &m_h;
    PRNG &m_prng;
    const size_t m_nspinorb;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
    const defs::inds m_exsigs;

public:
    ExcitGen(const Hamiltonian &h, PRNG &prng, defs::inds exsigs);

    virtual bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    virtual bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn);

    virtual size_t approx_nconn() const {
        return 0ul;
    }

    virtual std::string description() const {
        std::vector<std::string> strings;
        for (auto const &exsig: m_exsigs) strings.push_back(exsig_utils::to_string(exsig));
        return log::format("excitation generator for exsigs {}", string_utils::join(strings, ","));
    }

    virtual ~ExcitGen() {}
};

/**
 * Base class for stochastic Fermion sector excitations
 */
class FrmExcitGen : public ExcitGen {
protected:
    const bool m_spin_conserving;

public:
    FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t nexcit);

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override;
};

/**
 * Base class for stochastic Boson sector excitations
 */
//class BosExcitGen : public ExcitGen {
//protected:
//    const size_t m_nexcit;
//    const bool m_spin_conserving;
//
//public:
//    BosExcitGen(const Hamiltonian &h, PRNG &prng, size_t nexcit);
//
//    bool draw(const FrmOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
//              defs::ham_t &helem, conn::FrmOnv &conn) override;
//
//    bool draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
//              defs::ham_t &helem, conn::FrmBosOnv &conn) override;
//};

/**
 * Base class for stochastic excitations which may couple Fermion and Boson sectors excitations
 */
class FrmBosExcitGen : public ExcitGen {
public:
    FrmBosExcitGen(const Hamiltonian &h, PRNG &prng) :
            ExcitGen(h, prng, {exsig_utils::ex_1110, exsig_utils::ex_1101}) {}

    bool draw(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
              defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return ExcitGen::draw(exsig, src, orbs, prob, helem, conn);
    }

};

#endif //M7_EXCITGEN_H
