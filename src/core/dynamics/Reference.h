//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include "src/core/parallel/RankAllocator.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/Reducible.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "WalkerTable.h"
#include "Wavefunction.h"

class Reference : public elements::Onv {
    Wavefunction &m_wf;
    const Hamiltonian &m_ham;
    size_t m_irow;
    size_t m_irank;

    mutable conn::AsOnv m_aconn;
    Reducible<defs::ham_t> m_proj_energy_num;
    Reducible<defs::wf_comp_t> m_nwalker_at_doubles;
    Reducible<defs::wf_t> m_weight;

    /*
     * If a candidate for redefinition of the reference is found, then
     * its weight and row within m_list must be stored
     */
    Reducible<defs::wf_comp_t> m_candidate_weight;
    size_t m_irow_candidate;
    const double m_redefinition_thresh;
    bool m_redefinition_cycle;

public:
    Reference(Wavefunction &wf, const Hamiltonian& ham, views::Onv &onv, const Options &opts);

    void add_row(const size_t& irow){
        auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
        if (m_wf.m_walkers.m_flags.m_reference_connection(irow)) {
            add_to_numerator(m_wf.m_walkers.m_onv(irow), weight);
        }
        if (is_mine() && m_irow==irow) m_weight = weight;
    }

    using FermionOnv::operator=;
    void change(const size_t& irow, const size_t& irank){
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

    void log_candidate_weight(const size_t& irow, const defs::wf_comp_t& candidate_weight){
        if (irow==m_irow && mpi::i_am(m_irank)) return;
        if (candidate_weight>m_candidate_weight.local()){
            m_candidate_weight = candidate_weight;
            m_irow_candidate = irow;
        }
    }

    void reset() {
        m_proj_energy_num = 0.0;
        m_nwalker_at_doubles = 0.0;
        m_candidate_weight = 0.0;
        m_irow_candidate = ~0ul;
    }

    void reduce() {
        if (in_redefinition_cycle()) m_redefinition_cycle = false;
        m_candidate_weight.mpi_maxloc();
        m_weight.mpi_bcast(m_irank);
        if (m_candidate_weight.reduced()/std::abs(m_weight.reduced()) > m_redefinition_thresh) {
            //TODO: helpful output
//            std::cout<< " wwwww     " << m_candidate_weight.reduced() << "         " << std::abs(weight()) << std::endl;
//            std::cout<<m_candidate_weight.reduced()<< "  "<< m_candidate_weight.irank() << std::endl;
            if (mpi::i_am(m_candidate_weight.irank())){
//                std::cout << m_walkers.m_onv(m_irow_candidate).to_string() << std::endl;
            }
            change(m_irow_candidate, m_candidate_weight.irank());
        }
        m_proj_energy_num.mpi_sum();
        m_nwalker_at_doubles.mpi_sum(); // includes reference weight
    }

    const size_t &irow() {
        return m_irow;
    }

    const bool &in_redefinition_cycle() {
        return m_redefinition_cycle;
    }

    bool is_mine() const {
        return mpi::i_am(m_irank);
    }

    bool is_connected(const views::Onv &onv) const {
        auto &aconn = m_aconn;
        aconn.connect(*this, onv);
        return aconn.nexcit() < 3;
    }

    void add_to_numerator(const views::Onv &onv, const defs::wf_t &weight) {
        m_aconn.connect(*this, onv);
        m_proj_energy_num += m_ham.get_element(m_aconn) * weight;
        m_nwalker_at_doubles += std::abs(weight);
    }

    Reducible<defs::wf_comp_t> &nwalker_at_doubles() {
        return m_nwalker_at_doubles;
    }
    Reducible<defs::wf_comp_t> &candidate_weight() {
        return m_candidate_weight;
    }

    defs::ham_t proj_energy_num() const {
        return m_proj_energy_num.reduced();
    }
    defs::ham_comp_t weight() const {
        return m_weight.reduced();
    }
    defs::ham_comp_t proj_energy() const {
        return consts::real(proj_energy_num()/weight());
    }
};

#endif //M7_REFERENCE_H