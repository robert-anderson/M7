//
// Created by rja on 03/04/2022.
//

#ifndef M7_EXCITGEN2_H
#define M7_EXCITGEN2_H

#include <utility>

#include "M7_lib/sample/PRNG.h"
#include "M7_lib/field/Fields.h"
#include "M7_lib/connection/Connections.h"

struct ExcitGen2 {
    PRNG &m_prng;
    const defs::inds m_exsigs;
    const std::string m_description;
    ExcitGen2(PRNG &prng, defs::inds exsigs, std::string description):
        m_prng(prng), m_exsigs(std::move(exsigs)), m_description(std::move(description)){}

    virtual ~ExcitGen2() = default;

    /*
     * when the H matrix element is not necessary:
     */
    virtual bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

    virtual bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, conn::FrmBosOnv &conn);

    virtual bool draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn);

    /*
     * when the H matrix element is necessary. these can delegate the above methods in this base class, but in derived
     * classes it may make more sense to call specific methods to compute the matrix element in a more efficient way
     */
    virtual bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src,
                            defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) = 0;

    virtual bool draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) = 0;

    virtual bool draw_h_bos(const size_t &exsig, const field::BosOnv &src,
                            defs::prob_t &prob, defs::ham_t &helem, conn::BosOnv &conn) = 0;


    /*
     * here are defined homogeneously-named, statically defined dispatchers for the heterogeneously-named virtual
     * functions above. all draw_* functions could have been implemented as an overloaded virtual draw method if it were
     * not for the standard requirement that methods cannot be partially overridden. This dispatcher approach helps
     * cut down on clutter in the derived classes
     */
    bool draw(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
        return draw_frm(exsig, src, prob, conn);
    }
    bool draw(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, conn::FrmBosOnv &conn) {
        return draw_frmbos(exsig, src, prob, conn);
    }
    bool draw(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) {
        return draw_bos(exsig, src, prob, conn);
    }

    bool draw(const size_t &exsig, const field::FrmOnv &src,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    bool draw(const size_t &exsig, const field::FrmBosOnv &src,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn);

    bool draw(const size_t &exsig, const field::BosOnv &src,
              defs::prob_t &prob, defs::ham_t &helem, conn::BosOnv &conn);

    virtual size_t approx_nconn() const {
        return 1ul;
    }
};


#endif //M7_EXCITGEN2_H
