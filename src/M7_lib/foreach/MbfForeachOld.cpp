//
// Created by anderson on 2/9/22.
//

#include "MbfForeachOld.h"

mbf_foreach_old::frm::Base::Base(size_t nsite, body_fn_t body_fn, field_t *mbf) :
        mbf_foreach_old::Base<defs::Frm>({nsite, 0ul}, std::move(body_fn), mbf) {}

mbf_foreach_old::frm::Base::Base(const Base &other, field_t *mbf) :
        Base(other.m_bd.m_nsite, other.m_body_fn, mbf) {}

mbf_foreach_old::frm::General::Foreach::Foreach(mbf_foreach_old::frm::General &context, size_t nelec) :
        Ordered<>(context.m_mbf->m_nspinorb, nelec), m_context(context) {}


void mbf_foreach_old::frm::General::Foreach::body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) {
    m_context.m_mbf->zero();
    *m_context.m_mbf = value;
    m_context.body();
}

mbf_foreach_old::frm::General::General(size_t nsite, size_t nelec, body_fn_t body_fn, field_t *mbf) :
        Base(nsite, std::move(body_fn), mbf), m_foreach(*this, nelec) {}

mbf_foreach_old::frm::General::General(const General &other, field_t *mbf) :
        General(other.m_bd.m_nsite, other.m_foreach.m_nind, other.m_body_fn, mbf) {}

void mbf_foreach_old::frm::General::throwing_loop() {
    m_foreach.throwing_loop();
}

size_t mbf_foreach_old::frm::General::iiter() const {
    return m_foreach.iiter();
}

size_t mbf_foreach_old::frm::General::niter() const {
    return m_foreach.niter();
}

mbf_foreach_old::frm::Spins::Foreach::Foreach(mbf_foreach_old::frm::Spins &context, int ms2) :
        Ordered<>(context.m_mbf->m_nsite, (context.m_mbf->m_nsite + ms2) / 2), m_context(context) {}

void mbf_foreach_old::frm::Spins::Foreach::body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) {
    m_context.m_mbf->set_spins(value);
    m_context.body();
}

int mbf_foreach_old::frm::Spins::ms2() const {
    return 2 * m_foreach.m_nind - m_mbf->m_nsite;
}

mbf_foreach_old::frm::Spins::Spins(size_t nsite, int ms2, body_fn_t body_fn,
                                   field_t *mbf) :
        Base(nsite, std::move(body_fn), mbf), m_foreach(*this, ms2) {}

mbf_foreach_old::frm::Spins::Spins(const Spins &other, field_t *mbf) :
        Spins(other.m_bd.m_nsite, other.ms2(), other.m_body_fn, mbf) {}

void mbf_foreach_old::frm::Spins::throwing_loop() {
    m_foreach.throwing_loop();
}

size_t mbf_foreach_old::frm::Spins::iiter() const {
    return m_foreach.iiter();
}

size_t mbf_foreach_old::frm::Spins::niter() const {
    return m_foreach.niter();
}

mbf_foreach_old::frm::Ms2Conserve::Foreach::Foreach(mbf_foreach_old::frm::Ms2Conserve &context, size_t nelec) :
        Ordered<>(context.m_mbf->m_nsite, nelec), m_context(context) {}

void mbf_foreach_old::frm::Ms2Conserve::BetaForeach::body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) {
    m_context.m_mbf->put_spin_channel(1, false);
    m_context.m_mbf->set(m_context.m_mbf->m_nsite, value);
    m_context.body();
}

void mbf_foreach_old::frm::Ms2Conserve::AlphaForeach::body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) {
    m_context.m_mbf->put_spin_channel(0, false);
    m_context.m_mbf->set(0, value);
    m_context.m_beta_foreach.throwing_loop();
}

size_t mbf_foreach_old::frm::Ms2Conserve::nelec() const {
    return m_alpha_foreach.m_nind + m_beta_foreach.m_nind;
}

int mbf_foreach_old::frm::Ms2Conserve::ms2() const {
    return m_alpha_foreach.m_nind - m_beta_foreach.m_nind;
}

mbf_foreach_old::frm::Ms2Conserve::Ms2Conserve(size_t nsite, size_t nelec, int ms2, body_fn_t body_fn,
                                               field_t *mbf) :
        Base(nsite, std::move(body_fn), mbf),
        m_alpha_foreach(*this, nalpha(nelec, ms2)),
        m_beta_foreach(*this, nbeta(nelec, ms2)) {}

mbf_foreach_old::frm::Ms2Conserve::Ms2Conserve(const Ms2Conserve &other, field_t *mbf) :
        Ms2Conserve(other.m_bd.m_nsite, other.nelec(), other.ms2(), other.m_body_fn, mbf) {}

void mbf_foreach_old::frm::Ms2Conserve::throwing_loop() {
    m_alpha_foreach.throwing_loop();
}

size_t mbf_foreach_old::frm::Ms2Conserve::iiter() const {
    return m_alpha_foreach.iiter() * m_beta_foreach.niter() + m_beta_foreach.iiter();
}

size_t mbf_foreach_old::frm::Ms2Conserve::niter() const {
    return m_alpha_foreach.niter() * m_beta_foreach.niter();
}


mbf_foreach_old::bos::Base::Base(size_t nmode, body_fn_t body_fn, field_t *mbf) :
        mbf_foreach_old::Base<defs::Bos>({0ul, nmode}, std::move(body_fn), mbf) {}

mbf_foreach_old::bos::Base::Base(const Base &other, field_t *mbf) :
        Base(other.m_bd.m_nsite, other.m_body_fn, mbf) {}

mbf_foreach_old::bos::GeneralOpen::Foreach::Foreach(GeneralOpen &context, size_t nboson_max) :
        Unrestricted(context.m_mbf->m_nmode, nboson_max + 1), m_context(context) {}

void mbf_foreach_old::bos::GeneralOpen::Foreach::body(const inds_t &value, size_t iiter) {
    *m_context.m_mbf = value;
    m_context.body();
}

size_t mbf_foreach_old::bos::GeneralOpen::Foreach::nboson_max() const {
    if (m_shape.empty()) return ~0ul;
    return m_shape[0] - 1;
}


mbf_foreach_old::bos::GeneralOpen::GeneralOpen(size_t nmode, size_t nboson_max, body_fn_t body_fn, field::BosOnv *mbf) :
        Base(nmode, std::move(body_fn), mbf), m_foreach(*this, nboson_max) {}

mbf_foreach_old::bos::GeneralOpen::GeneralOpen(const mbf_foreach_old::bos::GeneralOpen &other, field::BosOnv *mbf) :
        GeneralOpen(other.m_bd.m_nmode, other.m_foreach.nboson_max(), other.m_body_fn, mbf) {}

void mbf_foreach_old::bos::GeneralOpen::throwing_loop() {
    m_foreach.throwing_loop();
}

size_t mbf_foreach_old::bos::GeneralOpen::iiter() const {
    return m_foreach.iiter();
}

size_t mbf_foreach_old::bos::GeneralOpen::niter() const {
    return m_foreach.niter();
}



void mbf_foreach_old::bos::GeneralClosed::throwing_loop() {
    m_foreach.throwing_loop();
}

size_t mbf_foreach_old::bos::GeneralClosed::iiter() const {
    return m_foreach.iiter();
}

size_t mbf_foreach_old::bos::GeneralClosed::niter() const {
    return m_foreach.niter();
}
