//
// Created by rja on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H


#include <src/core/sparse/SparseMatrix.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/parallel/Reducible.h>
#include <src/core/fermion/DeterminantList.h>
#include "WalkerList.h"

class DeterministicSubspace {

    WalkerList &m_walker_list;

    /*
     * indices of the deterministic subspace determinants in the
     * walker list on this process
     */
    struct SubspaceList : public DeterminantList {
        NumericField<size_t> m_irow; // row index in m_walker_list of DeterministicSubspace
        NumericField<size_t> m_irank;

        SubspaceList(std::string name, const size_t &nsite) : DeterminantList(name,1, nsite), m_irow(this), m_irank(this) {}
    };

    SubspaceList m_local_subspace_list;
    SubspaceList m_full_subspace_list;
    /*
     * just the weights stored on this MPI rank
     */
    std::vector<defs::wf_t> m_local_weights;
    /*
     * weights multiplied by H
     */
    std::vector<defs::wf_t> m_local_h_weights;
    /*
     * gathered weights from all processes
     */
    std::vector<defs::wf_t> m_full_weights;
    /*
     * sparse representation of H projected into the deterministic space
     */
    SparseMatrix<defs::ham_t> m_sparse_ham;

    defs::inds m_recvcounts;
    defs::inds m_displs;

    void build_hamiltonian(const Hamiltonian *ham) {
        ASSERT(m_sparse_ham.empty());
        m_full_subspace_list.all_gather(m_local_subspace_list);
        m_sparse_ham.resize(nrow_local());
        for (size_t irow_local = 0ul; irow_local < nrow_local(); ++irow_local) {
            // loop over local subspace (H rows)
            auto row_det = m_local_subspace_list.m_determinant(irow_local);
            for (size_t irow_full = 0ul; irow_full < nrow_full(); ++irow_full) {
                // loop over full subspace (H columns)
                // only add to sparse H if dets are connected
                auto col_det = m_full_subspace_list.m_determinant(irow_full);
                if (row_det == col_det) continue;
                auto helem = ham->get_element(row_det, col_det);
                if (!consts::float_is_zero(helem)) m_sparse_ham(irow_local, irow_full) = helem;
            }
        }

        m_local_weights.resize(nrow_local(), 0.0);
        m_local_h_weights.resize(nrow_local(), 0.0);
        m_full_weights.resize(nrow_full(), 0.0);
        mpi::all_gather(nrow_local(), m_recvcounts);
        mpi::counts_to_displs_consec(m_recvcounts, m_displs);
        ASSERT(m_displs.back() + m_recvcounts.back() == nrow_full());
    }

public:

    DeterministicSubspace(WalkerList &walker_list) :
            m_walker_list(walker_list),
            m_local_subspace_list(SubspaceList("local deterministic subspace", walker_list.m_determinant.m_nsite)),
            m_full_subspace_list(SubspaceList("full deterministic subspace", walker_list.m_determinant.m_nsite)),
            m_recvcounts(mpi::nrank(), 0), m_displs(mpi::nrank(), 0) {}

    void gather_and_project() {
        for (size_t irow_local = 0ul; irow_local < nrow_local(); ++irow_local) {
            auto irow_walker_list = *m_local_subspace_list.m_irow(irow_local);
            ASSERT(m_walker_list.m_flags.m_deterministic(irow_walker_list));
            m_local_weights[irow_local] = *m_walker_list.m_weight(irow_walker_list);
            ASSERT(m_local_weights[irow_local] == m_local_weights[irow_local])
        }
        ASSERT(nrow_local() == m_recvcounts[mpi::irank()])
        ASSERT(nrow_local() == m_local_weights.size())
        mpi::all_gatherv(m_local_weights, nrow_local(), m_full_weights, m_recvcounts, m_displs);
        m_local_h_weights.assign(m_local_h_weights.size(), 0);
        m_sparse_ham.multiply(m_full_weights, m_local_h_weights);
    }

    void rayleigh_quotient(defs::ham_t &num, defs::ham_comp_t &norm_square) {
        for (size_t irow_local = 0ul; irow_local < nrow_local(); ++irow_local) {
            auto w = m_local_weights[irow_local];
            auto hw = m_local_h_weights[irow_local];
            num += consts::conj(w) * hw;
            norm_square += std::pow(std::abs(w), 2.0);
        }
    }

    void update_weights(const double &tau, defs::wf_comp_t &delta_nw) {
        for (size_t irow_local = 0ul; irow_local < nrow_local(); ++irow_local) {
            auto irow_walker_list = *m_local_subspace_list.m_irow(irow_local);
            auto weight = m_walker_list.m_weight(irow_walker_list);
            ASSERT(*weight == *weight)
            auto h_weight = m_local_h_weights[irow_local];
            ASSERT(h_weight == h_weight)
            delta_nw -= std::abs(*weight);
            weight -= tau * h_weight;
            delta_nw += std::abs(*weight);
        }
    }

    void add_determinant(size_t irow_walker_list) {
        size_t irow = m_local_subspace_list.expand_push();
        m_local_subspace_list.m_determinant(irow) = m_walker_list.m_determinant(irow_walker_list);
        m_local_subspace_list.m_irow(irow) = irow_walker_list;
        m_walker_list.m_flags.m_deterministic(irow_walker_list) = true;
    }

    const size_t &nrow_local() const { return m_local_subspace_list.high_water_mark(0); }

    const size_t &nrow_full() const { return m_full_subspace_list.high_water_mark(0); }

    DeterminantElement local_det(const size_t &irow) const {
        ASSERT(irow < m_local_subspace_list.high_water_mark(0))
        return m_local_subspace_list.m_determinant(irow);
    }

    void build_from_whole_walker_list(Hamiltonian *ham) {
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(0); ++irow) {
            if (!m_walker_list.row_empty(irow)) add_determinant(irow);
        }
        build_hamiltonian(ham);
    }

    void build_from_det_connections(const DeterminantElement &ref, Hamiltonian *ham, size_t nexcit_max = 2) {
        std::cout << "Building deterministic subspace from connections of " << ref.to_string() << std::endl;
        Connection conn(ref);
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(0); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto det = m_walker_list.m_determinant(irow);
            conn.connect(ref, det);
            if (conn.nexcit() <= nexcit_max) add_determinant(irow);
        }
        build_hamiltonian(ham);
    }

    void build_from_nw_fraction(double fraction, defs::wf_comp_t nw, Hamiltonian *ham) {
        std::cout << "Building deterministic subspace of all dets with weight > total walker number * " << fraction << std::endl;
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(0); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto weight = m_walker_list.m_weight(irow);
            if (std::abs(*weight)/nw > fraction) add_determinant(irow);
        }
        build_hamiltonian(ham);
    }

    void build_from_highest_weighted(const WalkerList& list, size_t ndet_tot) {
        /*
        defs::inds row_inds;

        Connection conn(ref);
        for (size_t irow = 0ul; irow < m_walker_list.high_water_mark(0); ++irow) {
            if (m_walker_list.row_empty(irow)) continue;
            auto det = m_walker_list.m_determinant(irow);
            conn.connect(ref, det);
            if (conn.nexcit() <= nexcit_max) add_determinant(irow);
        }
        build_hamiltonian(ham);
         */
    }

};


#endif //M7_DETERMINISTICSUBSPACE_H
