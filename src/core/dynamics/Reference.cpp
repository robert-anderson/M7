//
// Created by rja on 03/07/2020.
//

#include "Reference.h"
#include "Wavefunction.h"

Reference::Reference(Wavefunction &wf, const Hamiltonian<> &ham, views::Onv<> &onv,
                     const Options &m_opts) :
        elements::Onv<>(onv), m_wf(wf), m_ham(ham), m_aconn(onv),
        m_redefinition_thresh(m_opts.reference_redefinition_thresh),
        m_proj_energy_num(m_summables, {1, 1}),
        m_nwalker_at_doubles(m_summables, {1, 1}),
        m_weight({1, 1}),
        m_candidate_weight({1, 1})
        {

    m_irank = m_wf.m_ra.get_rank(onv);
    /*
     * check for onv in walker table first
     */
    auto irow = wf.m_walkers[onv];
    m_irow = ~0ul;
    if (irow) {
        m_irow = *irow;
        ASSERT(m_irank == mpi::irank());
    } else if (mpi::i_am(m_irank)) {
        m_irow = m_wf.create_walker(onv, m_opts.nwalker_initial, ham.get_energy(onv), true);
    }
    change(m_irow, m_irank);
}

void Reference::add_row(const size_t &irow) {
    auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
    if (m_wf.m_walkers.m_flags.m_reference_connection(irow)) {
        add_to_numerator(m_wf.m_walkers.m_onv(irow), weight);
    }
    else {
        ASSERT(consts::float_is_zero(m_ham.get_element(*this, m_wf.m_walkers.m_onv(irow))))
    }
    if (is_mine() && m_irow==irow) m_weight(0, 0) = weight;
}

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
    if (candidate_weight>m_candidate_weight(0, 0)){
        m_candidate_weight(0, 0) = candidate_weight;
        m_irow_candidate = irow;
    }
}

void Reference::reset() {
    m_summables.zero();
}

void Reference::reduce() {
    if (in_redefinition_cycle()) m_redefinition_cycle = false;
    m_candidate_weight.all_maxloc();
    m_weight.bcast(m_irank);
    if (m_candidate_weight.reduced(0, 0)/std::abs(m_weight(0, 0)) > m_redefinition_thresh) {
        //TODO: helpful output
//            std::cout<< " wwwww     " << m_candidate_weight.reduced() << "         " << std::abs(weight()) << std::endl;
//            std::cout<<m_candidate_weight.reduced()<< "  "<< m_candidate_weight.irank() << std::endl;
        if (mpi::i_am(m_candidate_weight.rank_index(0, 0))){
//                std::cout << m_walkers.m_onv(m_irow_candidate).to_string() << std::endl;
        }
        change(m_irow_candidate, m_candidate_weight.rank_index(0, 0));
    }
    m_summables.all_sum();
}

const size_t &Reference::irow() {
    return m_irow;
}

const bool &Reference::in_redefinition_cycle() {
    return m_redefinition_cycle;
}

bool Reference::is_mine() const {
    return mpi::i_am(m_irank);
}

bool Reference::is_connected(const views::Onv<> &onv) const {
    m_aconn.connect(*this, onv);
    return m_aconn.connected();
}

void Reference::add_to_numerator(const views::Onv<> &onv, const defs::wf_t &weight) {
    m_aconn.connect(*this, onv);
    m_proj_energy_num(0, 0) += m_ham.get_element(m_aconn) * weight;
    m_nwalker_at_doubles(0, 0) += std::abs(weight);
}

ReductionMember<defs::wf_comp_t, defs::ndim_wf> &Reference::nwalker_at_doubles() {
    return m_nwalker_at_doubles;
}

ReductionMember<defs::wf_comp_t, defs::ndim_wf> &Reference::candidate_weight() {
    return m_candidate_weight;
}

defs::ham_t Reference::proj_energy_num() const {
    return m_proj_energy_num.reduced(0, 0);
}

defs::ham_comp_t Reference::weight() const {
    return m_weight(0, 0);
}

defs::ham_comp_t Reference::proj_energy() const {
    return consts::real(proj_energy_num()/weight());
}
