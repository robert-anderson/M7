//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITGEN_H
#define M7_EXCITGEN_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/connection/Connections.h>
#include "src/core/sample/PRNG.h"
#include "src/core/caches/CachedOrbs.h"

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
    virtual bool draw_frm(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::FrmOnv &conn);

    virtual bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::FrmBosOnv &conn);

    virtual bool draw_bos(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs,
                      defs::prob_t &prob, conn::BosOnv &conn);

    /*
     * when the H matrix element is necessary. these can delegate the above methods in this base class, but in derived
     * classes it may make more sense to call specific methods to compute the matrix element in a more efficient way
     */
    virtual bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                            defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    virtual bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
                               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn);

    virtual bool draw_h_bos(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs,
                            defs::prob_t &prob, defs::ham_t &helem, conn::BosOnv &conn);


    /*
     * here are defined homogeneously-named, statically defined dispatchers for the heterogeneously-named virtual
     * functions above. all draw_* functions could have been implemented as an overloaded virtual draw method if it were
     * not for the standard requirement that methods cannot be partially overridden. This dispatcher approach helps
     * cut down on clutter in the derived classes
     */
    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob, conn::FrmOnv &conn) {
        return draw_frm(exsig, src, orbs, prob, conn);
    }
    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob, conn::FrmBosOnv &conn) {
        return draw_frmbos(exsig, src, orbs, prob, conn);
    }
    bool draw(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob, conn::BosOnv &conn) {
        return draw_bos(exsig, src, orbs, prob, conn);
    }

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
        return draw_h_frm(exsig, src, orbs, prob, helem, conn);
    }
    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
        return draw_h_frmbos(exsig, src, orbs, prob, helem, conn);
    }
    bool draw(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::BosOnv &conn) {
        return draw_h_bos(exsig, src, orbs, prob, helem, conn);
    }

    virtual size_t approx_nconn() const {
        return 1ul;
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
    FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t exsig);

    bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, conn.m_frm);
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
