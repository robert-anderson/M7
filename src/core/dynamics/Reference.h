//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H


#include <src/core/parallel/RankAllocator.h>
#include <src/core/fermion/Determinant.h>
#include <src/core/fermion/Connection.h>
#include <src/core/parallel/Reducable.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include "WalkerList.h"

class Reference : public Determinant {

    WalkerList &m_list;
    RankAllocator<DeterminantElement> &m_ra;
    size_t m_irow;
    Reducable<size_t> m_irank;

    mutable AntisymConnection m_aconn;
    Reducable<defs::ham_t> m_proj_energy_num;
    Reducable<defs::wf_t> m_weight;

public:
    Reference(WalkerList &list, RankAllocator<DeterminantElement> &ra, DeterminantElement &det) :
            Determinant(det), m_list(list), m_ra(ra), m_aconn(det) {
        if (m_list.nrow_per_segment() == 0) m_list.expand(1);
        if (mpi::i_am(ra.get_rank(det))) {
            m_irow = m_list.push(det);
            ASSERT(m_irow != ~0ul)
            list.m_determinant(m_irow) = det;

            m_irank = mpi::irank();
        } else {
            m_irow = ~0ul;
            m_irank = ~0ul;
        }
        m_irank.mpi_min();
    }

    void update() {
        m_proj_energy_num = 0.0;
        if (is_mine()) m_weight = *m_list.m_weight(m_irow);
        m_weight.mpi_bcast(m_irank.reduced());
    }

    void synchronize() {
        m_proj_energy_num.mpi_sum();
    }

    const size_t &irow() {
        return m_irow;
    }

    bool is_mine() const {
        return mpi::i_am(m_irank.reduced());
    }

    bool is_connected(DeterminantElement &det) const {
        auto &aconn = m_aconn;
        aconn.connect(*this, det);
        return aconn.nexcit() < 3;
    }

    void add_to_numerator(const Hamiltonian *ham, const DeterminantElement &det, const defs::wf_t &weight) {
        m_aconn.connect(*this, det);
        m_proj_energy_num += ham->get_element(m_aconn) * weight;
    }

    Reducable<defs::ham_t> &proj_energy_num() {
        return m_proj_energy_num;
    }

    const defs::wf_t &weight() {
        return m_weight.reduced();
    }

    defs::ham_comp_t proj_energy() {
        return consts::real(proj_energy_num().reduced() / weight());
    }
};


#endif //M7_REFERENCE_H
