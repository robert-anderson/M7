//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"
#include "M7_lib/integrals/IntegralArray1e.h"
#include "M7_lib/integrals/IntegralArray2e.h"
#include "M7_lib/io/FcidumpFileReader.h"

struct GeneralFrmHam : FrmHam {

    typedef integrals_1e::Array<defs::ham_t> ints_1e_t;
    typedef integrals_2e::Array<defs::ham_t> ints_2e_t;
    typedef std::unique_ptr<ints_1e_t> ints_1e_ptr_t;
    typedef std::unique_ptr<ints_2e_t> ints_2e_ptr_t;

    struct Integrals {
        ints_1e_ptr_t m_1e;
        ints_2e_ptr_t m_2e;
    };

private:
    const FcidumpInfo m_info;

    static void log_ints_sym(integrals_1e::syms::Sym sym, bool initial);

    static void log_ints_sym(integrals_2e::syms::Sym sym, bool initial);

    Integrals make_ints(const FcidumpInfo& info, bool spin_major);

public:

    const Integrals m_ints;

    GeneralFrmHam(const FcidumpInfo& info, bool spin_major);

    explicit GeneralFrmHam(opt_pair_t opts);

    defs::ham_t get_coeff_1100(size_t a, size_t i) const override;

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

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

    size_t default_nelec() const override;

    int default_ms2_value() const override;
};


#endif //M7_GENERALFRMHAM_H
