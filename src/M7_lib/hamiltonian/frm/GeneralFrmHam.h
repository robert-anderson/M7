//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"

struct GeneralFrmHam : FrmHam {

    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    GeneralFrmHam(const sys::frm::Basis& hs);

    GeneralFrmHam(const FcidumpHeader& header, bool spin_major, int ms2_restrict=~0, size_t nelec=0);

    explicit GeneralFrmHam(const conf::FrmHam &opts);

    defs::ham_t get_coeff_1100(size_t a, size_t i) const override;

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    buffered::FrmOnv guess_reference() const;

    /**
     * @param prng
     *  random number generator
     * @param opts
     *  propagator options (so excitation generation settings can be read from the configuration document)
     * @param h
     *  Fermion H object to be **treated** as a GeneralFrmHam
     * @return
     *  list of excit gens required to stochastically propagate the given h
     */
    static excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& opts, const FrmHam& h);

    excit_gen_list_t make_excit_gens(PRNG &prng, const conf::Propagator& opts) const override;

    conn_foreach_list_t make_foreach_iters() const override;
};


#endif //M7_GENERALFRMHAM_H
