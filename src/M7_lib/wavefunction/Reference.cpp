//
// Created by Robert J. Anderson on 03/07/2020.
//

#include "Reference.h"
#include "Wavefunction.h"

wf::Ref::Ref(const conf::Reference &opts, const Hamiltonian &ham,
             const wf::Vectors &wf, uint_t ipart, TableBase::Loc loc) :
        shared_rows::Walker("reference", wf.m_store, loc),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_conn(ham.m_basis.size()),
        m_redefinition_thresh(opts.m_redef_thresh){
    if (m_redefinition_thresh==0.0)
        logging::info("Reference redefinition is deactivated");
    else
        REQUIRE_GE_ALL(m_redefinition_thresh, 1.0, "invalid redefinition threshold");
    logging::info("Initial reference MBF for WF part {} is {} with energy {}",
                  m_ipart, mbf(), gathered().m_row.m_hdiag);
}

void wf::Ref::update_ref_conn_flags() {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row; ++row){
        if (row.m_mbf.is_zero()) continue;
        row.m_ref_conn.put(m_ipart, is_connected(row.m_mbf));
    }
}

void wf::Ref::accept_candidate(uint_t icycle) {
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
    if (std::abs(candidate_weight) > std::abs(weight() * m_redefinition_thresh)){
        logging::info("Changing reference for WF part {} on cycle {}", m_ipart, icycle);
        logging::info("Current: {}, weight: {: .6e}, MPI rank: {}",
                      mbf().to_string(), weight(), m_wf.m_dist.irank(mbf()));
        redefine({irank_with_candidate, m_irow_candidate});
        logging::info("New    : {}, weight: {: .6e}, MPI rank: {}",
                      mbf().to_string(), candidate_weight, m_wf.m_dist.irank(mbf()));
        m_candidate_weight = 0.0;
        update_ref_conn_flags();
    }
}

void wf::Ref::contrib_row(const ::Walker& walker) {
    auto weight = 0.5*(walker.m_weight[m_ipart]+walker.m_weight[m_wf.ipart_replica(m_ipart)]);
    if (std::abs(weight) > std::abs(m_candidate_weight)) {
        m_candidate_weight = std::abs(weight);
        m_irow_candidate = walker.index();
    }
    if (walker.m_ref_conn.get(m_ipart)) {
        make_numerator_contribs(walker.m_mbf, walker.m_weight[m_ipart]);
    }
}

void wf::Ref::begin_cycle(uint_t icycle) {
    accept_candidate(icycle);
    m_candidate_weight = 0.0;
    m_proj_energy_num.m_local.zero();
    update();
}

void wf::Ref::end_cycle(uint_t /*icycle*/) {
    m_proj_energy_num.all_sum();
}

bool wf::Ref::is_connected(const field::Mbf &mbf) const {
    m_conn[mbf].connect(this->mbf(), mbf);
    return ham::is_significant(m_ham.get_element(this->mbf(), m_conn[mbf]));
}

void wf::Ref::make_numerator_contribs(const field::Mbf &mbf, const wf_t& weight) {
    m_conn[mbf].connect(mbf, this->mbf());
    m_proj_energy_num.m_local += m_ham.get_element(mbf, m_conn[mbf]) * weight;
}

const ham_t& wf::Ref::proj_energy_num() const {
    return m_proj_energy_num.m_reduced;
}

wf::Refs::Refs(const conf::Reference &opts, const Hamiltonian &ham, const wf::Vectors &wf, v_t<TableBase::Loc> locs) :
        m_proj_energy_nums(wf.m_format.m_shape), m_weights(wf.m_format.m_shape){
    DEBUG_ASSERT_EQ(locs.size(), wf.m_format.m_nelement,
                    "there should be a parallel table location specifying each reference row");
    m_refs.reserve(wf.m_format.m_nelement);
    for (uint_t ipart=0ul; ipart<wf.m_format.m_nelement; ++ipart) m_refs.emplace_back(opts, ham, wf, ipart, locs[ipart]);
}

const wf::Ref & wf::Refs::operator[](const uint_t &ipart) const {
    DEBUG_ASSERT_LT(ipart, m_refs.size(), "reference part index OOB");
    return m_refs[ipart];
}

void wf::Refs::begin_cycle(uint_t icycle) {
    for (auto& ref: m_refs) ref.begin_cycle(icycle);
}

void wf::Refs::end_cycle(uint_t icycle) {
    for (auto& ref: m_refs) ref.end_cycle(icycle);
}

void wf::Refs::contrib_row(const Walker& walker) {
    for (auto& ref: m_refs) ref.contrib_row(walker);
}

v_t<bool> wf::Refs::is_connected(const field::Mbf &mbf) const {
    v_t<bool> out;
    out.reserve(m_refs.size());
    for (uint_t ipart=0ul; ipart<m_refs.size(); ++ipart)
        out.push_back(m_refs[ipart].is_connected(mbf));
    return out;
}

const field::Numbers<ham_t, c_ndim_wf> & wf::Refs::proj_energy_nums() {
    uint_t ipart = 0ul;
    for (auto& ref: m_refs) m_proj_energy_nums[ipart++] = ref.proj_energy_num();
    return m_proj_energy_nums;
}

const field::Numbers<wf_t, c_ndim_wf> & wf::Refs::weights() {
    uint_t ipart = 0ul;
    for (auto& ref: m_refs) {
        m_weights[ipart] = ref.weight();
        ++ipart;
    }
    return m_weights;
}
