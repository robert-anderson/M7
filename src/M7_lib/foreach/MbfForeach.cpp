//
// Created by rja on 27/03/2022.
//

#include "MbfForeach.h"

mbf_foreach::Base::Base(const HilbertSpace& hs, size_t niter) :
    m_hs(hs), m_niter(niter), m_mbfs(m_hs) {}

void mbf_foreach::Base::loop(field::FrmOnv &mbf, const mbf_foreach::Base::function_t<field::FrmOnv> &fn) {
    frm_loop(mbf, fn);
}

void mbf_foreach::Base::loop(const mbf_foreach::Base::function_t<field::FrmOnv> &fn) {
    loop(m_mbfs.m_row.m_frm, fn);
}

void mbf_foreach::Base::loop(field::BosOnv& mbf, const mbf_foreach::Base::function_t<field::BosOnv> &fn) {
    bos_loop(mbf, fn);
}
void mbf_foreach::Base::loop(const mbf_foreach::Base::function_t<field::BosOnv> &fn) {
    loop(m_mbfs.m_row.m_bos, fn);
}

void mbf_foreach::Base::loop(field::FrmBosOnv& mbf, const mbf_foreach::Base::function_t<field::FrmBosOnv> &fn) {
    frmbos_loop(mbf, fn);
}
void mbf_foreach::Base::loop(const mbf_foreach::Base::function_t<field::FrmBosOnv> &fn) {
    loop(m_mbfs.m_row.m_frmbos, fn);
}

mbf_foreach::PairBase::PairBase(size_t nrow) : m_nrow(nrow), m_niter(nrow * nrow) {}

void mbf_foreach::PairBase::loop(const mbf_foreach::PairBase::function_t<field::FrmOnv> &fn) { frm_loop(fn); }

void mbf_foreach::PairBase::loop(const mbf_foreach::PairBase::function_t<field::BosOnv> &fn) { bos_loop(fn); }

void mbf_foreach::PairBase::loop(const mbf_foreach::PairBase::function_t<field::FrmBosOnv> &fn) { frmbos_loop(fn); }

mbf_foreach::frm::Base::Base(const FrmHilbertSpace& hs, size_t niter): mbf_foreach::Base({hs, {}}, niter) {}

mbf_foreach::frm::General::General(const FrmHilbertSpace& hs):
        Base(hs, foreach_t::niter(hs.m_sites.m_nspinorb, hs.m_nelec)),
        m_foreach(hs.m_sites.m_nspinorb, hs.m_nelec) {}

void mbf_foreach::frm::General::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm::Spins::Spins(const FrmHilbertSpace& hs) :
        Base(hs, foreach_t::niter(hs.m_sites, hs.m_nelec_alpha)),
        m_foreach(hs.m_sites, hs.m_nelec_alpha){}

void mbf_foreach::frm::Spins::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

size_t mbf_foreach::frm::Ms2Conserve::niter(const FrmHilbertSpace& hs) {
    auto na = foreach_t::niter(hs.m_sites, hs.m_nelec_alpha);
    auto nb = foreach_t::niter(hs.m_sites, hs.m_nelec_beta);
    return na * nb;
}

mbf_foreach::frm::Ms2Conserve::Ms2Conserve(const FrmHilbertSpace& hs) :
    Base(hs, niter(hs)),
    m_alpha_foreach(hs.m_sites, hs.m_nelec_alpha), m_beta_foreach(hs.m_sites, hs.m_nelec_beta){}

void mbf_foreach::frm::Ms2Conserve::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::Base::Base(const BosHilbertSpace& hs, size_t niter) :
    mbf_foreach::Base({{}, hs}, niter){}

mbf_foreach::bos::GeneralClosed::GeneralClosed(const BosHilbertSpace& hs) :
        Base(hs, foreach_t::niter(hs.m_nmode, hs.m_nboson)), m_foreach(hs.m_nmode, hs.m_nboson) {
    REQUIRE_TRUE(hs.m_nboson_conserve, "Not a closed quantum system in the bosonic sector");
}

void mbf_foreach::bos::GeneralClosed::bos_loop(field::BosOnv &mbf, const std::function<void(const field::BosOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::GeneralOpen::GeneralOpen(const BosHilbertSpace& hs):
        Base(hs, foreach_t::niter(hs.m_nmode, hs.m_occ_cutoff + 1)),
        m_foreach(hs.m_nmode, hs.m_occ_cutoff + 1) {}

void mbf_foreach::bos::GeneralOpen::bos_loop(field::BosOnv &mbf, const std::function<void(const field::BosOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm_bos::Base::Base(const HilbertSpace& hs, size_t niter) : mbf_foreach::Base(hs, niter) {}
