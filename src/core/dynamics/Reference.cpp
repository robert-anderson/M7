//
// Created by rja on 03/07/2020.
//

#include "Reference.h"

Reference::Reference(const Options &m_opts, const Hamiltonian<> &ham,
                     const Wavefunction& wf, size_t ipart, TableBaseZ::Loc loc):
        Wavefunction::DynamicRow(wf, loc, "reference"),
        m_ham(ham), m_wf(wf), m_ipart(ipart), m_aconn(ham.nsite()),
        m_redefinition_thresh(m_opts.reference_redefinition_thresh),
        m_proj_energy_num(m_summables, wf.m_format),
        m_nwalker_at_doubles(m_summables, wf.m_format)
        {
            update();
}

void Reference::add_row() {
    auto& row = m_wf.m_store.m_row;
    auto weight = row.m_weight(m_ipart);
    if (std::abs(weight) > m_candidate_abs_weight){
        m_candidate_abs_weight = std::abs(weight);
        m_irow_candidate = row.m_i;
    }
    if (row.m_reference_connection.get(m_ipart)) {
        add_to_numerator(row.m_onv, weight);
    }
}

#if 0
void Reference::change(const size_t &irow, const size_t &irank) {
    ASSERT(irank<mpi::nrank())
    ASSERT(irow!=~0ul || !mpi::i_am(irank))
    if (mpi::i_am(irank)){
        *this = m_wf.m_walkers.m_onv(irow);
    }
    /*
     * send this new reference definition to all MPI ranks
     */
    mpi_bcast(irank);
    /*
     * check that the storage of the new reference is consistent with dynamic
     * rank allocation
     */
    ASSERT(m_wf.m_ra.get_rank(*this)==irank)
    /*
     * the next cycle will re-evaluate the reference connection flag in the
     * Wavefunction object's m_walker_list member
     */
    m_redefinition_cycle = true;
    std::cout << "Reference ONV " << to_string() << " stored on MPI rank " << irank << std::endl;
    m_irank = irank;
    m_irow = irow;
}

void Reference::log_candidate_weight(const size_t &irow, const defs::wf_comp_t &candidate_weight) {
    if (irow==m_irow && mpi::i_am(m_irank)) return;
    if (candidate_weight>m_candidate_abs_weight(0, 0)){
        m_candidate_abs_weight(0, 0) = candidate_weight;
        m_irow_candidate = irow;
    }
}
#endif

void Reference::begin_cycle() {
    m_summables.zero();
    update();
}

void Reference::end_cycle() {
#if 0
    if (in_redefinition_cycle()) m_redefinition_cycle = false;
    m_candidate_abs_weight.all_maxloc();
    m_weight.bcast(m_irank);
    if (m_candidate_abs_weight.reduced(0, 0)/std::abs(m_weight(0, 0)) > m_redefinition_thresh) {
        if (mpi::i_am(m_candidate_abs_weight.rank_index(0, 0))){
        }
        change(m_irow_candidate, m_candidate_abs_weight.rank_index(0, 0));
    }
#endif
    m_summables.all_sum();
}

/*
const bool &Reference::in_redefinition_cycle() {
    return m_redefinition_cycle;
}*/

bool Reference::is_connected(const fieldsz::Onv<> &onv) const {
    m_aconn.connect(m_onv, onv);
    return m_aconn.connected();
}

void Reference::add_to_numerator(const fieldsz::Onv<> &onv, const defs::wf_t &weight) {
    m_aconn.connect(m_onv, onv);
    m_proj_energy_num(0, 0) += m_ham.get_element(m_aconn) * weight;
    m_nwalker_at_doubles(0, 0) += std::abs(weight);
}

ReductionMember<defs::wf_comp_t, defs::ndim_wf> &Reference::nwalker_at_doubles() {
    return m_nwalker_at_doubles;
}

defs::ham_t Reference::proj_energy_num() const {
    return m_proj_energy_num.reduced(0, 0);
}

defs::ham_comp_t Reference::proj_energy() const {
    return consts::real(proj_energy_num()/m_weight(m_ipart));
}

void Reference::update() {
    //accept_candidate(m_redefinition_thresh);
    Wavefunction::DynamicRow::update();
}
