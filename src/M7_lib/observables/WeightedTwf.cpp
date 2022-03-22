//
// Created by jhalson on 06/04/2021.
//

#include "WeightedTwf.h"

WeightedTwf::WeightedTwf(const Hamiltonian& ham, size_t npart, double_t fermion_factor, double_t boson_factor) :
        SpfTwfBase(ham, npart),
        m_frm_doub_occ_penalty_factor(fermion_factor),
        m_bos_occ_penalty_factor(boson_factor) {}

void WeightedTwf::add(const Numbers<defs::wf_t, defs::ndim_wf> &weight, defs::ham_t helem_sum, defs::ham_t diag_fac) {
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += weight[ipart] * helem_sum;
        m_denominator[ipart] += weight[ipart] * diag_fac;
    }
}

void WeightedTwf::add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                      const field::FrmOnv &onv) {
    auto diag_fac = evaluate_static_twf(onv);
    defs::ham_t helem_sum = diag_fac * m_ham.get_element(onv);
    auto fn = [&](const field::FrmOnv& dst, defs::ham_t helem){
        helem_sum+= evaluate_static_twf(dst)*helem;
    };
    m_excit_iters.foreach<field::FrmOnv>(onv, fn, true);
    add(weight, helem_sum, diag_fac);
}

void WeightedTwf::add(const Numbers<defs::wf_t, defs::ndim_wf> &weight, const field::FrmBosOnv &onv) {
    auto diag_fac = evaluate_static_twf(onv);
    defs::ham_t helem_sum = diag_fac * m_ham.get_element(onv);
    auto fn = [&](const field::FrmBosOnv& dst, defs::ham_t helem){
        helem_sum+= evaluate_static_twf(dst)*helem;
    };
    m_excit_iters.foreach<field::FrmBosOnv>(onv, fn, true);
    add(weight, helem_sum, diag_fac);
}

void WeightedTwf::add(const Numbers<defs::wf_t, defs::ndim_wf> &weight, const field::BosOnv &onv) {
    ABORT("not yet implemented");
}

defs::ham_t WeightedTwf::evaluate_static_twf(const field::FrmOnv &onv) const {
    size_t num_double_occ_sites = 0;
    for (size_t isite = 0; isite < onv.m_nsite; isite++) {
        num_double_occ_sites += (onv.get({0, isite}) and onv.get({1, isite}));
    }
    return std::exp(-m_frm_doub_occ_penalty_factor * num_double_occ_sites);
}

defs::ham_t WeightedTwf::evaluate_static_twf(const field::FrmBosOnv &onv) const {
    size_t total_boson_occ = 0;
    for (size_t isite = 0; isite < onv.m_frm.m_nsite; isite++) {
        total_boson_occ += onv.m_bos[isite];
    }
    return std::exp(-m_bos_occ_penalty_factor * total_boson_occ) + evaluate_static_twf(onv.m_frm);
}