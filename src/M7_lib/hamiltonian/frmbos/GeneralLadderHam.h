//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALLADDERHAM_H
#define M7_GENERALLADDERHAM_H

#include "M7_lib/io/EbdumpFileReader.h"
#include "M7_lib/integrals/FrmBosCoupledCoeffs.h"
#include "M7_lib/conf/HamiltonianConf.h"

#include "FrmBosHam.h"

struct GeneralLadderHam : FrmBosHam {
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

    GeneralLadderHam(const EbdumpHeader& header, const FrmHam& frm, const BosHam& bos, bool spin_major=false);

    GeneralLadderHam(const conf::FrmBosHam &opts, const FrmHam& frm, const BosHam& bos):
            GeneralLadderHam(EbdumpHeader(opts.m_ebdump.m_path), frm, bos, opts.m_ebdump.m_spin_major){}

    defs::ham_t get_coeff_0010(size_t imode) const override;

    defs::ham_t get_coeff_0001(size_t imode) const override;

    defs::ham_t get_coeff_1110(size_t imode, size_t i, size_t j) const override;

    defs::ham_t get_coeff_1101(size_t imode, size_t i, size_t j) const override;

    defs::ham_t get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const override;


    defs::ham_t get_element_pure(const field::FrmBosOnv &onv, size_t imode, bool cre) const;

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


    defs::ham_t get_element_coupled(const field::FrmBosOnv &onv,
                                    const conn::FrmOnv& frm_conn, size_t imode, bool cre) const;

    defs::ham_t get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


};

#endif //M7_GENERALLADDERHAM_H
