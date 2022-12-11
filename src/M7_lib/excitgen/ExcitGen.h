//
// Created by Robert J. Anderson on 03/04/2022.
//

#ifndef M7_EXCITGEN_H
#define M7_EXCITGEN_H

#include <utility>

#include "M7_lib/sample/PRNG.h"
#include "M7_lib/field/Fields.h"
#include "M7_lib/connection/Connections.h"

struct ExcitGen {
    PRNG& m_prng;
    const v_t<OpSig> m_exsigs;
    const str_t m_description;
    ExcitGen(PRNG& prng, v_t<OpSig> exsigs, str_t description):
        m_prng(prng), m_exsigs(std::move(exsigs)), m_description(std::move(description)){}

    virtual ~ExcitGen() = default;

    /*
     * when the H matrix element is not necessary:
     */
    virtual bool draw_frm(OpSig /*exsig*/, const field::FrmOnv& /*src*/, prob_t& prob, conn::FrmOnv& /*conn*/);

    virtual bool draw_frmbos(OpSig /*exsig*/, const field::FrmBosOnv& /*src*/, prob_t& prob, conn::FrmBosOnv& /*conn*/);

    virtual bool draw_bos(OpSig /*exsig*/, const field::BosOnv& /*src*/, prob_t& prob, conn::BosOnv& /*conn*/);

    /*
     * when the H matrix element is necessary. these can delegate the above methods in this base class, but in derived
     * classes it may make more sense to call specific methods to compute the matrix element in a more efficient way
     */
    virtual bool draw_h_frm(OpSig exsig, const field::FrmOnv& src,
                            prob_t& prob, ham_t& helem, conn::FrmOnv& conn) = 0;

    virtual bool draw_h_frmbos(OpSig exsig, const field::FrmBosOnv& src,
                               prob_t& prob, ham_t& helem, conn::FrmBosOnv& conn) = 0;

    virtual bool draw_h_bos(OpSig exsig, const field::BosOnv& src,
                            prob_t& prob, ham_t& helem, conn::BosOnv& conn) = 0;


    /*
     * get the probability of the excitgen drawing conn given the src MBF
     */
    virtual prob_t prob_frm(const field::FrmOnv& /*src*/, const conn::FrmOnv& /*conn*/) const {return 0.0;}
    virtual prob_t prob_bos(const field::BosOnv& /*src*/, const conn::BosOnv& /*conn*/) const {return 0.0;}
    virtual prob_t prob_frmbos(const field::FrmBosOnv& /*src*/, const conn::FrmBosOnv& /*conn*/) const {return 0.0;}

    /*
     * get the probability of the excitgen drawing conn given the src MBF and the ham matrix element already computed
     */
    virtual prob_t prob_h_frm(const field::FrmOnv& src, const conn::FrmOnv& conn, ham_t /*helem*/) const;
    virtual prob_t prob_h_bos(const field::BosOnv& src, const conn::BosOnv& conn, ham_t /*helem*/) const;
    virtual prob_t prob_h_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn, ham_t /*helem*/) const ;

    /*
     * here are defined homogeneously-named, statically defined dispatchers for the heterogeneously-named virtual
     * functions above. all draw_* functions could have been implemented as an overloaded virtual draw method if it were
     * not for the standard requirement that methods cannot be partially overridden. This dispatcher approach helps
     * cut down on clutter in the derived classes
     */
    bool draw(OpSig exsig, const field::FrmOnv& src, prob_t& prob, conn::FrmOnv& conn) {
        return draw_frm(exsig, src, prob, conn);
    }
    bool draw(OpSig exsig, const field::FrmBosOnv& src, prob_t& prob, conn::FrmBosOnv& conn) {
        return draw_frmbos(exsig, src, prob, conn);
    }
    bool draw(OpSig exsig, const field::BosOnv& src, prob_t& prob, conn::BosOnv& conn) {
        return draw_bos(exsig, src, prob, conn);
    }

    bool draw(OpSig exsig, const field::FrmOnv& src,
              prob_t& prob, ham_t& helem, conn::FrmOnv& conn);

    bool draw(OpSig exsig, const field::FrmBosOnv& src,
              prob_t& prob, ham_t& helem, conn::FrmBosOnv& conn);

    bool draw(OpSig exsig, const field::BosOnv& src,
              prob_t& prob, ham_t& helem, conn::BosOnv& conn);


    prob_t prob(const field::FrmOnv& src, const conn::FrmOnv& conn) const;
    prob_t prob(const field::BosOnv& src, const conn::BosOnv& conn) const;
    prob_t prob(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn) const;

    prob_t prob(const field::FrmOnv& src, const conn::FrmOnv& conn, ham_t helem) const;
    prob_t prob(const field::BosOnv& src, const conn::BosOnv& conn, ham_t helem) const;
    prob_t prob(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn, ham_t helem) const;

    virtual uint_t approx_nconn(OpSig /*exsig*/, sys::Particles /*particles*/) const {
        return 1ul;
    }

    typedef std::unique_ptr<ExcitGen> excit_gen_ptr_t;
    typedef std::forward_list<excit_gen_ptr_t> excit_gen_list_t;
};


#endif //M7_EXCITGEN_H
