//
// Created by rja on 03/07/2020.
//

#include "Reference.h"

Reference::Reference(const fciqmc_config::Reference &opts, const Hamiltonian &ham,
                     const Wavefunction &wf, size_t ipart, TableBase::Loc loc) :
        Wavefunction::SharedRow(wf, loc, "reference"),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_conn(ham.nsite()),
        m_redefinition_thresh(opts.m_redef_thresh){
    m_summables.add_members(m_proj_energy_num, m_nwalker_at_doubles);
    log::info("Initial reference ONV for WF part {} is {} with energy {}",
              m_ipart, get_mbf(), m_global.m_row.m_hdiag);
}

const field::Mbf &Reference::get_mbf() const {
    return m_global.m_row.m_mbf;
}

void Reference::update_ref_conn_flags() {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (row.m_mbf.is_zero()) continue;
        row.m_ref_conn.put(m_ipart, is_connected(row.m_mbf));
    }
}

void Reference::accept_candidate(double redefinition_thresh) {
    std::vector<defs::wf_comp_t> gather(mpi::nrank());
    mpi::all_gather(m_candidate_abs_weight, gather);
    DEBUG_ASSERT_EQ(m_candidate_abs_weight, gather[mpi::irank()], "Gather error");
    size_t irank = std::distance(gather.begin(), std::max_element(gather.begin(), gather.end()));
    mpi::bcast(m_irow_candidate, irank);
    auto current_weight = weight();
    if (std::abs(gather[irank]) > std::abs(current_weight*redefinition_thresh)){
        log::info("Changing the reference ONV for WF part {}. current ONV: {}, weight: {}",
                  m_ipart, get_mbf().to_string(), current_weight);
        redefine({irank, m_irow_candidate});
        log::info("Changed the reference ONV for WF part {}. new ONV: {}, weight: {}",
                  m_ipart, get_mbf().to_string(), gather[irank]);
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
        make_numerator_contribs(row.m_mbf, row.m_weight[m_ipart]);
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

bool Reference::is_connected(const field::Mbf &mbf) const {
    m_conn[mbf].connect(get_mbf(), mbf);
    return !consts::float_is_zero(m_ham.get_element(get_mbf(), m_conn[mbf]));
}

void Reference::make_numerator_contribs(const field::Mbf &mbf, const defs::wf_t& weight) {
    m_conn[mbf].connect(get_mbf(), mbf);
    m_proj_energy_num.m_local += m_ham.get_element(get_mbf(), m_conn[mbf]) * weight;
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
    auto unnorm = m_global.m_row.m_average_weight[ipart]+m_global.m_row.m_weight[ipart];
    return unnorm/static_cast<defs::wf_comp_t>(m_global.m_row.occupied_ncycle(icycle));
}

References::References(const fciqmc_config::Reference &opts, const Hamiltonian &ham, const Wavefunction &wf,
                       std::vector<TableBase::Loc> locs) :
        m_proj_energy_nums(wf.m_format.m_shape), m_weights(wf.m_format.m_shape){
    ASSERT(locs.size()==wf.m_format.m_nelement);
    m_refs.reserve(wf.m_format.m_nelement);
    for (size_t ipart=0ul; ipart<wf.m_format.m_nelement; ++ipart) m_refs.emplace_back(opts, ham, wf, ipart, locs[ipart]);
    ASSERT(m_refs.size()==wf.npart());
}

const Reference &References::operator[](const size_t &ipart) const {
    return m_refs[ipart];
}

void References::begin_cycle() {
    for (auto& ref: m_refs) ref.begin_cycle();
}

void References::end_cycle() {
    for (auto& ref: m_refs) ref.end_cycle();
}

void References::contrib_row() {
    for (auto& ref: m_refs) ref.contrib_row();
}

std::vector<bool> References::is_connected(const field::Mbf &mbf) const {
    std::vector<bool> out;
    out.reserve(m_refs.size());
    for (size_t ipart=0ul; ipart<m_refs.size(); ++ipart)
        out.push_back(m_refs[ipart].is_connected(mbf));
    return out;
}

const field::Numbers<defs::ham_t, defs::ndim_wf> &References::proj_energy_nums() {
    size_t ipart = 0ul;
    for (auto& ref: m_refs) m_proj_energy_nums[ipart++] = ref.proj_energy_num();
    return m_proj_energy_nums;
}

const field::Numbers<defs::wf_t, defs::ndim_wf> &References::weights() {
    size_t ipart = 0ul;
    for (auto& ref: m_refs) m_weights[ipart++] = ref.weight();
    return m_weights;
}
