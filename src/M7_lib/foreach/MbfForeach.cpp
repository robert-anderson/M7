//
// Created by Robert J. Anderson on 27/03/2022.
//

#include "MbfForeach.h"

mbf_foreach::Base::Base(size_t niter) : m_niter(niter) {}

mbf_foreach::PairBase::PairBase(size_t nrow) : m_nrow(nrow), m_niter(nrow * nrow) {}

mbf_foreach::frm::Base::Base(const sys::frm::Sector& sector, size_t niter):
    mbf_foreach::Base(niter), m_sector(sector) {}

mbf_foreach::frm::General::General(const sys::frm::Sector& sector):
        Base(sector, foreach_t::niter(sector.m_basis.m_nspinorb, sector.m_elecs)),
        m_foreach(sector.m_basis.m_nspinorb, sector.m_elecs) {}

void mbf_foreach::frm::General::frm_loop(field::FrmOnv& mbf, const function_t &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm::Spins::Spins(const sys::frm::Sector& sector) :
        Base(sector, foreach_t::niter(sector.m_basis.m_nsite, sector.m_elecs.m_nalpha)),
        m_foreach(sector.m_basis.m_nsite, sector.m_elecs.m_nalpha){}

void mbf_foreach::frm::Spins::frm_loop(field::FrmOnv& mbf, const function_t &fn) {
    loop_fn(mbf, fn);
}

size_t mbf_foreach::frm::Ms2Conserve::niter(const sys::frm::Sector& sector) {
    return sector.size();
}

mbf_foreach::frm::Ms2Conserve::Ms2Conserve(const sys::frm::Sector& sector) :
    Base(sector, niter(sector)),
    m_alpha_foreach(sector.m_basis.m_nsite, sector.m_elecs.m_nalpha),
    m_beta_foreach(sector.m_basis.m_nsite, sector.m_elecs.m_nbeta){
    REQUIRE_TRUE(sector.m_elecs.m_ms2.conserve(), "2*Ms must be conserved");
}

void mbf_foreach::frm::Ms2Conserve::frm_loop(field::FrmOnv& mbf, const function_t &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::Base::Base(const sys::bos::Sector& sector, size_t niter) :
    mbf_foreach::Base(niter), m_sector(sector){}

mbf_foreach::bos::GeneralClosed::GeneralClosed(const sys::bos::Sector& sector) :
        Base(sector, foreach_t::niter(sector.m_basis.m_nmode, sector.m_bosons)),
        m_foreach(sector.m_basis.m_nmode, sector.m_bosons) {
    REQUIRE_TRUE(sector.m_bosons.conserve(), "Not a closed quantum system in the bosonic sector");
}

void mbf_foreach::bos::GeneralClosed::bos_loop(field::BosOnv &mbf, const function_t &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::bos::GeneralOpen::GeneralOpen(const sys::bos::Sector& sector):
        Base(sector, foreach_t::niter(sector.m_basis.m_nmode, sector.m_basis.m_occ_cutoff + 1)),
        m_foreach(sector.m_basis.m_nmode, sector.m_basis.m_occ_cutoff + 1) {}

void mbf_foreach::bos::GeneralOpen::bos_loop(field::BosOnv &mbf, const function_t &fn) {
    loop_fn(mbf, fn);
}

mbf_foreach::frm_bos::Base::Base(const sys::Sector& sector, size_t niter) : mbf_foreach::Base(niter) {}
