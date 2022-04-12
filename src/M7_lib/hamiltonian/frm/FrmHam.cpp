//
// Created by rja on 27/02/2020.
//

#include "M7_lib/enumerator/ContainerCombinationEnumerator.h"
#include "M7_lib/caches/DecodedDeterminants.h"
#include "M7_lib/parallel/RankAllocator.h"

#include "FrmHam.h"

FrmHam::FrmHam(size_t nelec, size_t nsite, int ms2_restrict, const defs::inds& site_irreps):
        m_nelec(nelec), m_nsite(nsite), m_ms2_restrict(ms2_restrict),
        m_point_group_map(PointGroup(), site_irreps.empty() ? defs::inds(nsite, 0ul) : site_irreps),
        m_contribs_1100(exsig_utils::ex_single), m_contribs_2200(exsig_utils::ex_double) {
}

defs::ham_t FrmHam::get_element(const field::FrmOnv &onv) const {
    return get_element_0000(onv);
}

defs::ham_comp_t FrmHam::get_energy(const field::FrmOnv &onv) const {
    return consts::real(get_element_0000(onv));
}

defs::ham_t FrmHam::get_element(const field::FrmOnv &ket, const conn::FrmOnv &conn) const {
    switch (conn.size()) {
        case 0:
            return get_element_0000(ket);
        case 2:
            return get_element_1100(ket, conn);
        case 4:
            return get_element_2200(ket, conn);
        default:
            return 0.0;
    }
}

size_t FrmHam::nci() const {
    return ci_utils::fermion_dim(m_nsite, m_nelec);
}

void FrmHam::log_data() const {
    if (disabled()) return;
    if (!m_contribs_1100.is_nonzero(0ul))
        log::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(exsig_utils::ex_single))
        log::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(0ul))
        log::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_single))
        log::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_double))
        log::info("2-electron term has no double-excitation contributions");
}
