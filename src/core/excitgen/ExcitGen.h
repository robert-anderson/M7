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
    const BasisDims m_bd;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
public:
    const defs::inds m_exsigs;

    ExcitGen(const Hamiltonian &h, PRNG &prng, defs::inds exsigs);


    /*
     * when the H matrix element is not necessary:
     */
    virtual bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::FrmOnv &conn);

    virtual bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::FrmBosOnv &conn);

    virtual bool draw(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::BosOnv &conn);

    /*
     * when the H matrix element is necessary. these can delegate the above methods in this base class, but in derived
     * classes it may make more sense to call specific methods to compute the matrix element in a more efficient way
     */
    virtual bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    virtual bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn);

    virtual bool draw(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::BosOnv &conn);


    virtual size_t approx_nconn() const {
        return 0ul;
    }

    virtual std::string description() const = 0;

    virtual ~ExcitGen() {}
};

/**
 * Base class for stochastic Fermion sector excitations
 */
struct FrmExcitGen : public ExcitGen {
    using ExcitGen::draw;
protected:
    const bool m_spin_conserving;
public:
    FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t nexcit);


    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, conn.m_frm);
    }

    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, helem, conn.m_frm);
    }

};

/**
 * Base class for stochastic excitations which involve single boson creation or annihilation operators
 */
struct LadderExcitGen : public ExcitGen {
    LadderExcitGen(const Hamiltonian &h, PRNG &prng, defs::inds exsigs) :
            ExcitGen(h, prng, exsigs) {}
};


/**
 * Base class for stochastic Boson number-conserving excitations
 */
struct BosExcitGen : public ExcitGen {
    BosExcitGen(const Hamiltonian &h, PRNG &prng, size_t exsig);

    size_t approx_nconn() const override;
};

#endif //M7_EXCITGEN_H
