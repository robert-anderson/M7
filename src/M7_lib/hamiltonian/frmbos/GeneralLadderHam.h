//
// Created by Robert J. Anderson on 12/9/21.
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

    GeneralLadderHam(const EbdumpInfo& info, bool spin_major, uint_t bos_occ_cutoff=sys::bos::c_max_occ);

    GeneralLadderHam(opt_pair_t opts):
        GeneralLadderHam(EbdumpInfo(opts.m_ham.m_ebdump.m_path),
                opts.m_ham.m_ebdump.m_spin_major, opts.m_basis.m_bos_occ_cutoff){}

    defs::ham_t get_coeff_1110(uint_t imode, uint_t i, uint_t j) const override;

    defs::ham_t get_coeff_1101(uint_t imode, uint_t i, uint_t j) const override;


    defs::ham_t get_element_pure(const field::FrmBosOnv &onv, uint_t imode, bool cre) const;

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


    defs::ham_t get_element_coupled(const field::FrmBosOnv &onv,
                                    const conn::FrmOnv& frm_conn, uint_t imode, bool cre) const;

    defs::ham_t get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;


};

#endif //M7_GENERALLADDERHAM_H
