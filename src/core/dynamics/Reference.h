//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H


#include "src/core/parallel/RankAllocator.h"
#include "src/core/basis/Determinant.h"
#include "src/core/basis/Connection.h"
#include "src/core/parallel/Reducible.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "WalkerList.h"

class Reference : public Determinant {

    WalkerList &m_list;
    RankAllocator<DeterminantElement> &m_ra;
    size_t m_irow;
    size_t m_irank;

    mutable AntisymConnection m_aconn;
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
    Reference(WalkerList &list, RankAllocator<DeterminantElement> &ra,
            DeterminantElement &det, double redefinition_thresh) :
            Determinant(det), m_list(list), m_ra(ra), m_aconn(det),
            m_redefinition_thresh(redefinition_thresh){
        if (m_list.nrow_per_segment() == 0) m_list.expand(1);
        m_irank = ra.get_rank(det);
        if (mpi::i_am(m_irank)) {
            m_irow = m_list.push(det);
            ASSERT(m_irow != ~0ul)
            list.m_determinant(m_irow) = det;
        } else {
            m_irow = ~0ul;
        }
        change(m_irow, m_irank);
    }

    using Determinant::operator=;
    void change(const size_t& irow, const size_t& irank){
        ASSERT(irank<mpi::nrank())
        ASSERT(irow!=~0ul || !mpi::i_am(irank))
        if (mpi::i_am(irank)){
            *this = m_list.m_determinant(irow);
        }
        /*
         * send this new reference definition to all MPI ranks
         */
        mpi_bcast(irank);
        /*
         * check that the storage of the new reference is consistent with dynamic
         * rank allocation
         */
        ASSERT(m_ra.get_rank(*this)==irank)
        /*
         * the next cycle will re-evaluate the reference connection flag in the
         * Wavefunction object's m_walker_list member
         */
        m_redefinition_cycle = true;
        std::cout << "Reference determinant " << to_string() << " stored on MPI rank " << irank << std::endl;
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

    void update() {
        m_proj_energy_num = 0.0;
        m_nwalker_at_doubles = 0.0;
        m_candidate_weight = 0.0;
        m_irow_candidate = ~0ul;
        if (is_mine()) m_weight = *m_list.m_weight(m_irow);
        m_weight.mpi_bcast(m_irank);
    }

    void synchronize() {
        if (in_redefinition_cycle()) m_redefinition_cycle = false;
        m_candidate_weight.mpi_maxloc();
        if (m_candidate_weight.reduced()/std::abs(weight()) > m_redefinition_thresh) {
            //TODO: helpful output
            std::cout<< " wwwww     " << m_candidate_weight.reduced() << "         " << std::abs(weight()) << std::endl;
            std::cout<<m_candidate_weight.reduced()<< "  "<< m_candidate_weight.irank() << std::endl;
            if (mpi::i_am(m_candidate_weight.irank())){
                std::cout << m_list.m_determinant(m_irow_candidate).to_string() << std::endl;
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

    bool is_connected(const DeterminantElement &det) const {
        auto &aconn = m_aconn;
        aconn.connect(*this, det);
        return aconn.nexcit() < 3;
    }

    void add_to_numerator(const Hamiltonian *ham, const DeterminantElement &det, const defs::wf_t &weight) {
        m_aconn.connect(*this, det);
        m_proj_energy_num += ham->get_element(m_aconn) * weight;
        m_nwalker_at_doubles += std::abs(weight);
    }

    Reducible<defs::ham_t> &proj_energy_num() {
        return m_proj_energy_num;
    }
    Reducible<defs::wf_comp_t> &nwalker_at_doubles() {
        return m_nwalker_at_doubles;
    }
    Reducible<defs::wf_comp_t> &candidate_weight() {
        return m_candidate_weight;
    }

    const defs::wf_t &weight() {
        return m_weight.reduced();
    }

    defs::ham_comp_t proj_energy() {
        return consts::real(proj_energy_num().reduced() / weight());
    }
};


#endif //M7_REFERENCE_H
