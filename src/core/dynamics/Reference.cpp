//
// Created by rja on 03/07/2020.
//

#include "Reference.h"

Reference::Reference(const fciqmc_config::Reference &opts, const Hamiltonian<> &ham,
                     const Wavefunction &wf, size_t ipart, TableBase::Loc loc) :
        Wavefunction::SharedRow(wf, loc, "reference"),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_aconn(ham.nsite()),
        m_redefinition_thresh(opts.m_redef_thresh){
    m_summables.add_members(m_proj_energy_num, m_nwalker_at_doubles);
    log::debug("Initial reference ONV for WF part {} is {}", m_ipart, get_onv().to_string());
    ASSERT(wf.m_store.is_protected());
}

const fields::Onv<> &Reference::get_onv() const {
    return m_global.m_row.m_onv;
}

void Reference::update_ref_conn_flags() {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (row.m_onv.is_zero()) continue;
        row.m_ref_conn.put(m_ipart, is_connected(row.m_onv));
    }
}

void Reference::accept_candidate(double redefinition_thresh) {
    std::vector<defs::wf_t> gather(mpi::nrank());
    mpi::all_gather(m_candidate_abs_weight, gather);
    DEBUG_ASSERT_EQ(m_candidate_abs_weight, gather[mpi::irank()], "Gather error");
    size_t irank = std::distance(gather.begin(), std::max_element(gather.begin(), gather.end()));
    mpi::bcast(m_irow_candidate, irank);
    auto current_weight = weight();
    if (std::abs(gather[irank]) > std::abs(current_weight*redefinition_thresh)){
        log::debug("Changing the reference ONV for WF part {}. current ONV: {}, weight: {}",
                   m_ipart, get_onv().to_string(), current_weight);
        redefine({irank, m_irow_candidate});
        log::debug("Changed the reference ONV for WF part {}. new ONV: {}, weight: {}",
                   m_ipart, get_onv().to_string(), gather[irank]);
        m_candidate_abs_weight = 0.0;
        update_ref_conn_flags();
    }
}

void Reference::contrib_row() {
    auto &row = m_wf.m_store.m_row;
    auto weight = row.m_weight[m_ipart];
    if (std::abs(weight) > m_candidate_abs_weight) {
        m_candidate_abs_weight = std::abs(weight);
        m_irow_candidate = row.index();
    }
    if (row.m_ref_conn.get(m_ipart)) {
        make_numerator_contribs(row.m_onv, row.m_weight[m_ipart]);
    }
}

void Reference::begin_cycle() {
    m_candidate_abs_weight = 0.0;
    m_summables.zero_all_local();
    update();
}

void Reference::end_cycle() {
    accept_candidate(m_redefinition_thresh);
    m_summables.all_sum();
}

bool Reference::is_connected(const fields::Onv<> &onv) const {
    m_aconn.connect(get_onv(), onv);
    return !consts::float_is_zero(m_ham.get_element(m_aconn));
}

bool Reference::connection_phase(const fields::Onv<> &onv) const {
    m_aconn.connect(get_onv(), onv);
    return m_aconn.phase();
}

void Reference::make_numerator_contribs(const fields::Onv<> &onv, const defs::wf_t& weight) {
    m_aconn.connect(get_onv(), onv);
    m_proj_energy_num.m_local += m_ham.get_element(m_aconn) * weight;
    m_nwalker_at_doubles.m_local += std::abs(weight);
}

const defs::wf_comp_t& Reference::nwalker_at_doubles() {
    return m_nwalker_at_doubles.m_reduced;
}

const defs::ham_t& Reference::proj_energy_num() const {
    return m_proj_energy_num.m_reduced;
}


const defs::wf_t &Reference::weight() const {
    return m_global.m_row.m_weight[m_ipart];
}

defs::wf_t Reference::norm_average_weight(const size_t& icycle, const size_t& ipart) const {
    return (m_global.m_row.m_average_weight[ipart]+m_global.m_row.m_weight[ipart])/(m_global.m_row.occupied_ncycle(icycle));
}