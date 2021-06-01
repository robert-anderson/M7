//
// Created by rja on 03/07/2020.
//

#include "Reference.h"

Reference::Reference(const Options &m_opts, const Hamiltonian<> &ham,
                     const Wavefunction &wf, size_t ipart, TableBase::Loc loc) :
        Wavefunction::SharedRow(wf, loc, "reference"),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_aconn(ham.nsite()),
        m_redefinition_thresh(m_opts.reference_redefinition_thresh),
        m_proj_energy_num(wf.m_format),
        m_nwalker_at_doubles(wf.m_format) {
    m_summables.add_members(m_proj_energy_num, m_nwalker_at_doubles);
    ASSERT(wf.m_store.is_protected());
}

const fields::Onv<> &Reference::get_onv() const {
    return m_global.m_row.m_onv;
}

void Reference::update_ref_conn_flags() {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (row.m_onv.is_zero()) continue;
        row.m_reference_connection.put(m_ipart, is_connected(row.m_onv));
    }
}

void Reference::accept_candidate(double redefinition_thresh) {
    std::vector<defs::wf_t> gather(mpi::nrank());
    mpi::all_gather(m_candidate_abs_weight, gather);
    MPI_ASSERT(m_candidate_abs_weight==gather[mpi::irank()], "Gather error");
    size_t irank = std::distance(gather.begin(), std::max_element(gather.begin(), gather.end()));
    mpi::bcast(m_irow_candidate, irank);
    auto current_weight = weight()[m_ipart];
    if (gather[irank] > std::abs(current_weight*redefinition_thresh)){
        log::debug("Changing the reference ONV to {} for WF part {}. current weight: {}, "
                   "candidate weight: {}", get_onv().to_string(), m_ipart, current_weight, gather[irank]);
        redefine({irank, m_irow_candidate});
        m_candidate_abs_weight = 0.0;
        update_ref_conn_flags();
    }
}

void Reference::contrib_row() {
    auto &row = m_wf.m_store.m_row;
    auto weight = row.m_weight[m_ipart];
    if (std::abs(weight) > m_candidate_abs_weight) {
        m_candidate_abs_weight = std::abs(weight);
        m_irow_candidate = row.m_i;
    }
    if (row.m_reference_connection.get(m_ipart)) {
        make_numerator_contribs(row.m_onv, row.m_weight);
    }
}

void Reference::begin_cycle() {
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

void Reference::make_numerator_contribs(const fields::Onv<> &onv, const fields::Numbers<defs::ham_t, defs::ndim_wf>& weights) {
    m_aconn.connect(get_onv(), onv);
    m_proj_energy_num.m_local.add_scaled(m_ham.get_element(m_aconn), weights);
    m_nwalker_at_doubles.m_local.add_abs(weights);
}

NdReduction<defs::wf_comp_t, defs::ndim_wf> &Reference::nwalker_at_doubles() {
    return m_nwalker_at_doubles;
}

const fields::Numbers<defs::ham_t, defs::ndim_wf>& Reference::proj_energy_num() const {
    return m_proj_energy_num.m_reduced;
}


const fields::Numbers<defs::ham_t, defs::ndim_wf> &Reference::weight() const {
    return m_global.m_row.m_weight;
}

defs::wf_t Reference::norm_average_weight(const size_t& icycle, const size_t& ipart) const {
    return (m_global.m_row.m_average_weight[ipart]+m_global.m_row.m_weight[ipart])/(m_global.m_row.occupied_ncycle(icycle));
}