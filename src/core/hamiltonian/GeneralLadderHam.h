//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALLADDERHAM_H
#define M7_GENERALLADDERHAM_H

#include <io/EbdumpFileReader.h>
#include <integrals/FrmBosCoupledCoeffs.h>
#include <config/Hamiltonian.h>
#include "LadderHam.h"

struct GeneralLadderHam : LadderHam {
    /**
     * coefficients for "coupled" ranksigs 1110, 1101. contributing exsigs are either:
     *  "density coupled" (0010, 0001), or
     *  "hopping coupled" (1110, 1101)
     */
    FrmBosCoupledCoeffs m_v;
    /**
     * coefficients for "uncoupled" ranksigs 0010, 0001. only contributing exsigs are
     *  "uncoupled" (0010, 0001)
     *
     * density-coupled and uncoupled excitations have the same exsig, collectively they will be called "pure" bosonic
     * excitations / de-excitations, as opposed to the fermion-coupled "hopping" exsigs
     */
    std::vector<defs::ham_t> m_v_unc;

    GeneralLadderHam(const EbdumpHeader& header, size_t nboson_max);

    GeneralLadderHam(const fciqmc_config::LadderHamiltonian &opts) :
        GeneralLadderHam(EbdumpHeader(opts.m_ebdump.m_path), opts.m_nboson_max){}

    defs::ham_t get_coeff_0010(const size_t &imode) const override;

    defs::ham_t get_coeff_0001(const size_t &imode) const override;

    defs::ham_t get_coeff_1110(const size_t &imode, const size_t &j, const size_t &i) const override;

    defs::ham_t get_coeff_1101(const size_t &imode, const size_t &j, const size_t &i) const override;

    defs::ham_t get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const override;


    defs::ham_t get_element_pure(const field::FrmBosOnv &onv, const size_t& imode, bool cre) const;

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


    defs::ham_t get_element_coupled(const field::FrmBosOnv &onv,
                                    const conn::FrmOnv& frm_conn, const size_t& imode, bool cre) const;

    defs::ham_t get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


};

#endif //M7_GENERALLADDERHAM_H