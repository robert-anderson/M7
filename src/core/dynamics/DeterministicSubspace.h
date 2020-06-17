//
// Created by rja on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H


#include <src/core/sparse/SparseMatrix.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include "WalkerList.h"

class DeterministicSubspace {

    WalkerList* m_walker_list = nullptr;
    /*
     * indices of the deterministic subspace determinants in the
     * walker list on this process
     */
    struct SubspaceList: public List {
        DeterminantField m_determinant;
        NumericField<size_t> m_irow;
        NumericField<size_t> m_irank;
        SubspaceList(const size_t& nsite):
        m_determinant(this, 1, nsite), m_irow(this), m_irank(this) {}
    };
    std::unique_ptr<SubspaceList> m_local_subspace_list;
    std::unique_ptr<SubspaceList> m_full_subspace_list;
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

    void gather_weights(){
#pragma omp parallel for default(none)
        for (size_t irow_local=0ul; irow_local<nrow_local(); ++irow_local){
            auto irow_walker_list = *m_local_subspace_list->m_irow(irow_local);
            m_local_weights[irow_local] = *m_walker_list->m_weight(irow_walker_list);
        }
        mpi::all_gatherv(m_local_weights.data(), nrow_local(), m_full_weights.data(), m_recvcounts, m_displs);
    }

    void project_propagator(){
        m_local_h_weights.assign(m_local_h_weights.size(), 0);
        m_sparse_ham.multiply(m_local_weights, m_local_h_weights);
    }

public:
    DeterministicSubspace(WalkerList* walker_list):
    m_local_subspace_list(std::unique_ptr<SubspaceList>(new SubspaceList(walker_list->m_determinant.m_nsite))),
    m_full_subspace_list(std::unique_ptr<SubspaceList>(new SubspaceList(walker_list->m_determinant.m_nsite))),
    m_recvcounts(mpi::nrank(), 0), m_displs(mpi::nrank(), 0), m_walker_list(walker_list){}

    void add_determinant(size_t irow_walker_list){
        size_t irow = m_local_subspace_list->expand_push();
        m_local_subspace_list->m_determinant(irow) = m_walker_list->m_determinant(irow_walker_list);
    }

    const size_t& nrow_local() const {return m_local_subspace_list->high_water_mark(0);}
    const size_t& nrow_full() const {return m_full_subspace_list->high_water_mark(0);}

    void build_hamiltonian(const Hamiltonian* ham){
        ASSERT(m_sparse_ham.empty());
        m_full_subspace_list->all_gather(*m_local_subspace_list);
        m_sparse_ham.resize(nrow_local());
#pragma omp parallel for default(none) shared(ham)
        for (size_t irow_local=0ul; irow_local<nrow_local(); ++irow_local){
            // loop over local subspace (H rows)
            auto row_det = m_local_subspace_list->m_determinant(irow_local);
            for (size_t irow_full=0ul; irow_full<nrow_full(); ++irow_full){
                // loop over full subspace (H columns)
                // only add to sparse H if dets are connected
                auto col_det = m_full_subspace_list->m_determinant(irow_full);
                if (row_det==col_det) continue;
                auto helem = ham->get_element(row_det, col_det);
                if (!consts::float_is_zero(helem)) m_sparse_ham(irow_local, irow_full) = helem;
            }
        }
        m_local_weights.resize(nrow_local(), 0.0);
        m_local_h_weights.resize(nrow_local(), 0.0);
        m_full_weights.resize(nrow_full(), 0.0);
        size_t nrow = nrow_local();
        mpi::all_gather(&nrow, 1, m_recvcounts.data(), 1);
        for (size_t i=1ul; i<mpi::nrank(); ++i) m_displs[i] = m_displs[i-1]+m_recvcounts[i-1];
        ASSERT(m_displs.back()+m_recvcounts.back()==nrow_full());
    }

};


#endif //M7_DETERMINISTICSUBSPACE_H
