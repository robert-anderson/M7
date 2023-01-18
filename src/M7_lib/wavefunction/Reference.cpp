//
// Created by Robert J. Anderson on 03/07/2020.
//

#include "Reference.h"

Reference::Reference(const conf::Reference &opts, const Hamiltonian &ham,
                     const Wavefunction &wf, uint_t ipart, TableBase::Loc loc) :
        shared_rows::Walker("reference", wf.m_store, loc),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_conn(ham.m_basis.size()),
        m_redefinition_thresh(opts.m_redef_thresh){
    if (m_redefinition_thresh==0.0)
        logging::info("Reference redefinition is deactivated");
    else
        REQUIRE_GE_ALL(m_redefinition_thresh, 1.0, "invalid redefinition threshold");
    m_summables.add_members(m_proj_energy_num, m_nwalker_at_doubles);
    logging::info("Initial reference MBF for WF part {} is {} with energy {}",
              m_ipart, mbf(), m_all.m_row.m_hdiag);
}

void Reference::update_ref_conn_flags() {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row; ++row){
        if (row.m_mbf.is_zero()) continue;
        row.m_ref_conn.put(m_ipart, is_connected(row.m_mbf));
    }
}

void Reference::accept_candidate(uint_t icycle) {
    /*
     * only continue if the redefinition threshold is non-zero
     */
    if (m_redefinition_thresh==0.0) return;
    v_t<wf_comp_t> gather(mpi::nrank());
    mpi::all_gather(m_candidate_weight, gather);
    DEBUG_ASSERT_EQ(m_candidate_weight, gather[mpi::irank()], "Gather error");
    auto cmp_fn = [](const wf_t& w1, const wf_t& w2){return std::abs(w1) > std::abs(w2);};
    auto it_best = std::max_element(gather.cbegin(), gather.cend(), cmp_fn);
    DEBUG_ASSERT_FALSE(it_best==gather.cend(), "max element not found");
    uint_t irank_with_candidate = std::distance(gather.cbegin(), it_best);
    mpi::bcast(m_irow_candidate, irank_with_candidate);
    auto candidate_weight = *it_best;
    auto current_weight = weight(m_ipart);
    if (std::abs(candidate_weight) > std::abs(current_weight * m_redefinition_thresh)){
        logging::info("Changing reference for WF part {} on cycle {}", m_ipart, icycle);
        logging::info("Current: {}, weight: {: .6e}, MPI rank: {}",
                      mbf().to_string(), current_weight, m_wf.m_dist.irank(mbf()));
        redefine({irank_with_candidate, m_irow_candidate});
        logging::info("New    : {}, weight: {: .6e}, MPI rank: {}",
                      mbf().to_string(), candidate_weight, m_wf.m_dist.irank(mbf()));
        m_candidate_weight = 0.0;
        update_ref_conn_flags();
    }
}

void Reference::contrib_row(const ::Walker& walker) {
    auto weight = 0.5*(walker.m_weight[m_ipart]+walker.m_weight[m_wf.ipart_replica(m_ipart)]);
    if (std::abs(weight) > std::abs(m_candidate_weight)) {
        m_candidate_weight = std::abs(weight);
        m_irow_candidate = walker.index();
    }
    if (walker.m_ref_conn.get(m_ipart)) {
        make_numerator_contribs(walker.m_mbf, walker.m_weight[m_ipart]);
    }
}

void Reference::begin_cycle(uint_t icycle) {
    accept_candidate(icycle);
    m_candidate_weight = 0.0;
    m_summables.zero_all_local();
    update();
}

void Reference::end_cycle(uint_t /*icycle*/) {
    m_summables.all_sum();
}

bool Reference::is_connected(const field::Mbf &mbf) const {
    m_conn[mbf].connect(this->mbf(), mbf);
    return ham::is_significant(m_ham.get_element(this->mbf(), m_conn[mbf]));
}

void Reference::make_numerator_contribs(const field::Mbf &mbf, const wf_t& weight) {
    m_conn[mbf].connect(mbf, this->mbf());
    m_proj_energy_num.m_local += m_ham.get_element(mbf, m_conn[mbf]) * weight;
    m_nwalker_at_doubles.m_local += std::abs(weight);
}

const wf_comp_t& Reference::nwalker_at_doubles() {
    return m_nwalker_at_doubles.m_reduced;
}

const ham_t& Reference::proj_energy_num() const {
    return m_proj_energy_num.m_reduced;
}

References::References(const conf::Reference &opts, const Hamiltonian &ham, const Wavefunction &wf,
                       v_t<TableBase::Loc> locs) :
        m_proj_energy_nums(wf.m_format.m_shape), m_weights(wf.m_format.m_shape){
    DEBUG_ASSERT_EQ(locs.size(), wf.m_format.m_nelement,
                    "there should be a parallel table location specifying each reference row");
    m_refs.reserve(wf.m_format.m_nelement);
    for (uint_t ipart=0ul; ipart<wf.m_format.m_nelement; ++ipart) m_refs.emplace_back(opts, ham, wf, ipart, locs[ipart]);
}

const Reference &References::operator[](const uint_t &ipart) const {
    DEBUG_ASSERT_LT(ipart, m_refs.size(), "reference part index OOB");
    return m_refs[ipart];
}

void References::begin_cycle(uint_t icycle) {
    for (auto& ref: m_refs) ref.begin_cycle(icycle);
}

void References::end_cycle(uint_t icycle) {
    for (auto& ref: m_refs) ref.end_cycle(icycle);
}

void References::contrib_row(const Walker& walker) {
    for (auto& ref: m_refs) ref.contrib_row(walker);
}

v_t<bool> References::is_connected(const field::Mbf &mbf) const {
    v_t<bool> out;
    out.reserve(m_refs.size());
    for (uint_t ipart=0ul; ipart<m_refs.size(); ++ipart)
        out.push_back(m_refs[ipart].is_connected(mbf));
    return out;
}

const field::Numbers<ham_t, c_ndim_wf> &References::proj_energy_nums() {
    uint_t ipart = 0ul;
    for (auto& ref: m_refs) m_proj_energy_nums[ipart++] = ref.proj_energy_num();
    return m_proj_energy_nums;
}

const field::Numbers<wf_t, c_ndim_wf> &References::weights() {
    uint_t ipart = 0ul;
    for (auto& ref: m_refs) {
        m_weights[ipart] = ref.weight(ipart);
        ++ipart;
    }
    return m_weights;
}
