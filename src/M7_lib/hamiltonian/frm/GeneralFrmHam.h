//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"
#include "M7_lib/integrals/IntegralArray1e.h"
#include "M7_lib/integrals/IntegralArray2e.h"
#include "M7_lib/integrals/IntegralReader.h"
#include "M7_lib/io/FcidumpTextFileReader.h"

struct GeneralFrmHam : FrmHam {

    typedef integrals_1e::Array<ham_t> ints_1e_t;
    typedef integrals_2e::Array<ham_t> ints_2e_t;
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

    Integrals make_ints(IntegralReader& reader);

    Integrals make_ints(const FcidumpInfo& info, bool spin_major) {
        const str_t fmt = "Reading fermion Hamiltonian coefficients from {} file \"" + info.m_fname + "\"...";
        if (info.m_impl==FcidumpInfo::CSV) {
            log::info(fmt, "plain text CSV");
            CsvIntegralReader reader(info, spin_major);
            return make_ints(reader);
        }
        else if (info.m_impl==FcidumpInfo::MolcasHDF5) {
            log::info(fmt, "Molcas HDF5 binary");
            MolcasHdf5IntegralReader reader(info, spin_major);
            return make_ints(reader);
        }
        return {nullptr, nullptr};
    }

public:

    const Integrals m_ints;

    GeneralFrmHam(const FcidumpInfo& info, bool spin_major);

    explicit GeneralFrmHam(opt_pair_t opts);

    ham_t get_coeff_1100(uint_t a, uint_t i) const override;

    ham_t get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const override;


    ham_t get_element_0000(const field::FrmOnv &onv) const override;

    ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

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

    uint_t default_nelec() const override;

    int default_ms2_value() const override;
};


#endif //M7_GENERALFRMHAM_H
