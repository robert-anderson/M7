//
// Created by rja on 27/03/2022.
//

#include "MbfForeach.h"

mbf_foreach::Base::Base(BasisData bd, size_t niter) :
    m_bd(std::move(bd)), m_niter(niter), m_mbfs(m_bd) {}

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

mbf_foreach::frm::Base::Base(size_t nsite, size_t niter) : mbf_foreach::Base({nsite, 0ul}, niter) {}

mbf_foreach::frm::NumberConserve::NumberConserve(size_t nsite, size_t nelec, size_t niter) : Base(nsite, niter), m_nelec(nelec) {}

mbf_foreach::frm::General::General(size_t nsite, size_t nelec) :
        NumberConserve(nsite, nelec, foreach_t::niter(2 * nsite, nelec)),
        m_foreach(2 * nsite, nelec) {}

void mbf_foreach::frm::General::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm::Spins::Spins(size_t nsite, int ms2) :
        NumberConserve(nsite, nsite + ms2 / 2, foreach_t::niter(nsite, (nsite + ms2) / 2)),
        m_ms2(ms2), m_foreach(nsite, (nsite + ms2) / 2) {}

void mbf_foreach::frm::Spins::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

size_t mbf_foreach::frm::Ms2Conserve::niter(size_t nsite, size_t nelec, int ms2) {
    auto na = foreach_t::niter(nsite, ci_utils::nalpha(nelec, ms2));
    auto nb = foreach_t::niter(nsite, ci_utils::nbeta(nelec, ms2));
    return na * nb;
}

mbf_foreach::frm::Ms2Conserve::Ms2Conserve(size_t nsite, size_t nelec, int ms2) :
        NumberConserve(nsite, nelec, niter(nsite, nelec, ms2)), m_ms2(ms2),
        m_alpha_foreach(m_bd.m_nsite, ci_utils::nalpha(m_nelec, m_ms2)),
        m_beta_foreach(m_bd.m_nsite, ci_utils::nbeta(m_nelec, m_ms2)) {}

void mbf_foreach::frm::Ms2Conserve::frm_loop(field::FrmOnv& mbf, const std::function<void(const field::FrmOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::Base::Base(size_t nmode, size_t niter) : mbf_foreach::Base({0ul, nmode}, niter){}

mbf_foreach::bos::NumberConserve::NumberConserve(size_t nmode, size_t nboson, size_t niter) : Base(nmode, niter), m_nboson(nboson) {}

mbf_foreach::bos::GeneralClosed::GeneralClosed(size_t nmode, size_t nboson) :
        NumberConserve(nmode, nboson, foreach_t::niter(nmode, nboson)), m_foreach(nmode, nboson) {}

void mbf_foreach::bos::GeneralClosed::bos_loop(field::BosOnv &mbf, const std::function<void(const field::BosOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::GeneralOpen::GeneralOpen(size_t nmode, size_t nboson_max) :
        Base(nmode, foreach_t::niter(nmode, nboson_max + 1)),
        m_foreach(nmode, nboson_max + 1) {}

void mbf_foreach::bos::GeneralOpen::bos_loop(field::BosOnv &mbf, const std::function<void(const field::BosOnv &)> &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm_bos::Base::Base(BasisData bd, size_t niter) : mbf_foreach::Base(bd, niter) {}
