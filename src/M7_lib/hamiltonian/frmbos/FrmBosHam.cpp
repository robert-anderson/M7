//
// Created by rja on 05/11/2020.
//

#include "FrmBosHam.h"

#include <utility>

FrmBosHam::FrmBosHam(const HilbertSpace &hs, const FrmHam &frm, const BosHam &bos) :
        m_hs(FrmHilbertSpace(hs.m_frm, frm.m_hs), BosHilbertSpace(hs.m_bos, bos.m_hs)),
        m_contribs_0010(exsig_utils::ex_0010), m_contribs_0001(exsig_utils::ex_0001),
        m_contribs_1110(exsig_utils::ex_1110), m_contribs_1101(exsig_utils::ex_1101) {}

FrmBosHam::FrmBosHam(const FrmHam &frm, const BosHam &bos) : FrmBosHam({frm.m_hs, bos.m_hs}, frm, bos){}

void FrmBosHam::log_data() const {
    if (disabled()) return;
    if (!m_contribs_0010.is_nonzero(exsig_utils::ex_0010))
        log::info("0010 uncoupled boson ladder hamiltonian term has no contributions");
    if (!m_contribs_0001.is_nonzero(exsig_utils::ex_0001))
        log::info("0001 uncoupled boson ladder hamiltonian term has no contributions");

    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_0010))
        log::info("1110 fermion-coupled boson ladder term has no 0010 contributions");
    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_1110))
        log::info("1110 fermion-coupled boson ladder term has no 1110 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_0001))
        log::info("1101 fermion-coupled boson ladder term has no 0001 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_1101))
        log::info("1101 fermion-coupled boson ladder term has no 1101 contributions");
}