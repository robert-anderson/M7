//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H


#include <src/core/parallel/RankAllocator.h>
#include <src/core/fermion/Determinant.h>
#include <src/core/fermion/Connection.h>
#include <src/core/thread/PrivateStore.h>
#include <src/core/parallel/Hybrid.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include "WalkerList.h"

class Reference : public Determinant {

    WalkerList &m_list;
    RankAllocator<DeterminantElement> &m_ra;
    size_t m_irow;
    Distributed<size_t> m_irank;

    std::unique_ptr<PrivateStore<AntisymConnection>> m_aconn;
    Hybrid<defs::ham_t> m_proj_energy_num;
    Distributed<defs::wf_t> m_weight;

public:
    Reference(WalkerList &list, RankAllocator<DeterminantElement> &ra, DeterminantElement &det) :
            Determinant(det), m_list(list), m_ra(ra) {
        if (mpi::i_am(ra.get_rank(det))) {
            m_irow = list.push(det);
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
        m_proj_energy_num.put_thread_sum();
        m_proj_energy_num.mpi_sum();
    }

    const size_t &irow() {
        return m_irow;
    }

    bool is_mine() const {
        return mpi::i_am(m_irank.reduced());
    }

    bool is_connected(DeterminantElement &det) const {
        auto &aconn = m_aconn->get();
        aconn.connect(*this, det);
        return aconn.nexcit() < 3;
    }

    void add_to_numerator(const Hamiltonian *ham, const DeterminantElement &det, const defs::wf_t &weight) {
        auto &aconn = m_aconn->get();
        aconn.connect(*this, det);
        m_proj_energy_num.thread() += ham->get_element(aconn) * weight;
    }

    Hybrid<defs::ham_t>& proj_energy_num() {
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
