//
// Created by rja on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H

#include "src/core/sparse/SparseMatrix.h"
#include "src/core/hamiltonian/FermionHamiltonian.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/Reducible.h"
#include "src/core/basis/DeterminantList.h"
#include "WalkerTable.h"

class DeterministicSubspace {

    WalkerTable &m_walker_list;

    /*
     * indices of the deterministic subspace determinants in the
     * walker list on this process
     */
    struct SubspaceRow : Row {
        fields::FermionOnv m_onv;
        fields::Number<size_t> m_irow_wf; // row index in m_walker_list of DeterministicSubspace
        fields::Number<size_t> m_irank;

        SubspaceRow(const size_t &nsite) :
                Row(), m_onv(this, nsite), m_irow_wf(this), m_irank(this) {}
    };

    typedef Table<SubspaceRow> subspace_list_t;

    typedef std::vector<defs::wf_t> weight_vec_t;
    typedef std::vector<weight_vec_t> weight_vecs_t;

    /**
     * indices of WF parts subject to deterministic projection in this subspace
     */
    defs::inds m_subject_iparts;

    subspace_list_t m_local_subspace_list;
    subspace_list_t m_global_subspace_list;
    
    /*
     * just the weights stored on this MPI rank
     */
    weight_vecs_t m_local_weights;
    /*
     * weights multiplied by H
     */
    weight_vecs_t m_local_h_weights;
    /*
     * gathered weights from all processes
     */
    weight_vecs_t m_global_weights;
    /*
     * sparse representation of H projected into the deterministic space
     */
    sparse::Matrix<defs::ham_t> m_sparse_ham;

    defs::inds m_recvcounts;
    defs::inds m_displs;

    size_t nsubject() const {
        return m_subject_iparts.size();
    }

    void build_hamiltonian(const FermionHamiltonian &ham) {
        ASSERT(m_sparse_ham.empty());
        m_global_subspace_list.all_gatherv(m_local_subspace_list);
        std::cout << "Building semistochastic sparse hamiltonian with " << nrow_global() << " total rows and "
            << nrow_local() << " local rows" << std::endl;
        m_sparse_ham.resize(nrow_local());
        auto& row_local = m_local_subspace_list.m_row;
        for (row_local.restart(); row_local.in_range(); row_local.step()){
            // loop over local subspace (H rows)
            auto& row_global = m_global_subspace_list.m_row;
            for (row_global.restart(); row_global.in_range(); row_global.step()){
                // loop over full subspace (H columns)
                // only add to sparse H if dets are connected
                if (row_local.m_onv == row_global.m_onv) continue;
                auto helem = ham.get_element(row_local.m_onv, row_global.m_onv);
                if (!consts::float_is_zero(helem)) m_sparse_ham.add(row_local.m_i, row_global.m_i, helem);
            }
        }

        for (auto& v: m_local_weights) v.resize(nrow_local(), 0.0);
        for (auto& v: m_local_h_weights) v.resize(nrow_local(), 0.0);
        for (auto& v: m_global_weights) v.resize(nrow_global(), 0.0);
        mpi::all_gather(nrow_local(), m_recvcounts);
        mpi::counts_to_displs_consec(m_recvcounts, m_displs);
        ASSERT(m_displs.back() + m_recvcounts.back() == nrow_global());
    }

public:

    DeterministicSubspace(WalkerTable &walker_list, defs::inds subject_iparts) :
            m_walker_list(walker_list), m_subject_iparts(subject_iparts),
            m_local_subspace_list(SubspaceRow(walker_list.m_row.m_onv.m_nsite)),
            m_global_subspace_list(SubspaceRow(walker_list.m_row.m_onv.m_nsite)),
            m_local_weights(nsubject()), m_local_h_weights(nsubject()), m_global_weights(nsubject()),
            m_recvcounts(mpi::nrank(), 0), m_displs(mpi::nrank(), 0) {}

    void gather_and_project() {
        auto row_wf = m_walker_list.m_row;
        auto& row_local = m_local_subspace_list.m_row;
        for (row_local.restart(); row_local.in_range(); row_local.step()){
            row_wf.jump(row_local.m_irow_wf);
            ASSERT(row_wf.m_deterministic.get(0));
            for (size_t isubject=0ul; isubject<nsubject(); ++isubject)
                m_local_weights[row_local.m_i][isubject] = row_wf.m_weight[m_subject_iparts[isubject]];
        }
        ASSERT(nrow_local() == m_recvcounts[mpi::irank()])
        ASSERT(nrow_local() == m_local_weights.size())

        for (size_t isubject=0ul; isubject<nsubject(); ++isubject) {
            mpi::all_gatherv(m_local_weights[isubject], nrow_local(),
                             m_global_weights[isubject], m_recvcounts, m_displs);
            m_local_h_weights[isubject].assign(m_local_h_weights[isubject].size(), 0);
            m_sparse_ham.multiply(m_global_weights[isubject], m_local_h_weights[isubject]);
        }
    }

//    void rayleigh_quotient(defs::ham_t &num, defs::ham_comp_t &norm_square) {
//        for (size_t irow_local = 0ul; irow_local < nrow_local(); ++irow_local) {
//            auto w = m_local_weights[irow_local];
//            auto hw = m_local_h_weights[irow_local];
//            num += consts::conj(w) * hw;
//            norm_square += std::pow(std::abs(w), 2.0);
//        }
//    }

    void update_weights(const double &tau, fields::Numbers<defs::wf_t, defs::ndim_wf>& delta_nw) {
        auto row_wf = m_walker_list.m_row;
        auto& row_local = m_local_subspace_list.m_row;
        for (row_local.restart(); row_local.in_range(); row_local.step()){
            row_wf.jump(row_local.m_irow_wf);
            for (size_t isubject=0ul; isubject<nsubject(); ++isubject) {
                const auto& subject = m_subject_iparts[isubject];
                auto& weight = row_wf.m_weight[subject];
                const auto& h_weight = m_local_h_weights[isubject][row_local.m_i];
                delta_nw[subject] -= std::abs(weight);
                weight -= tau * h_weight;
                delta_nw[subject] += std::abs(weight);
            }
        }
    }

    void add_wf_row(WalkerTableRow& row_wf) {
        size_t irow = m_local_subspace_list.push_back();
        auto& row = m_local_subspace_list.m_row;

        row.jump(irow);
        row.m_onv = row_wf.m_onv;
        row.m_irow_wf = row_wf.m_i;

        for (size_t isubject=0ul; isubject<nsubject(); ++isubject) {
            const auto &subject = m_subject_iparts[isubject];
            row_wf.m_deterministic.set(subject);
        }
    }

    const size_t &nrow_local() const { return m_local_subspace_list.m_hwm;}

    const size_t &nrow_global() const { return m_global_subspace_list.m_hwm;}


    void build_from_whole_walker_list(const FermionHamiltonian &ham) {
        auto row = m_walker_list.m_row;
        for (row.restart(); row.in_range(); row.step()){
            if (!row.is_cleared()) add_wf_row(row);
        }
        build_hamiltonian(ham);
    }

    void build_from_onv_connections(const fields::Onv<> &ref, const FermionHamiltonian &ham, size_t nexcit_max) {
        std::cout << "Building deterministic subspace from connections of " << ref.to_string() << std::endl;
        conn::Antisym<> conn(ref);
        auto row = m_walker_list.m_row;
        for (row.restart(); row.in_range(); row.step()){
            if (row.is_cleared()) continue;
            conn.connect(ref, row.m_onv);
            if (conn.nexcit() <= nexcit_max) add_wf_row(row);
        }
        build_hamiltonian(ham);
    }

    void build_from_nw_fraction(double fraction, defs::wf_comp_t nw, const FermionHamiltonian &ham) {
        log::info("Building deterministic subspace of all dets with "
                  "weight > total walker number * {} in all WF parts", fraction);
        auto row = m_walker_list.m_row;
        auto large_enough = [&](const size_t& isubject){return std::abs(row.m_weight[isubject])/nw > fraction;};
        for (row.restart(); row.in_range(); row.step()){
            if (row.is_cleared()) continue;
            if (std::all_of(m_subject_iparts.cbegin(), m_subject_iparts.cend(), large_enough)) add_wf_row(row);
        }
        build_hamiltonian(ham);
    }


    void build_from_nadd_thresh(double thresh, defs::wf_comp_t nadd, const FermionHamiltonian &ham) {
        log::info("Building deterministic subspace of all dets with "
                  "weight > nadd ({}) * {} in all WF parts", nadd, thresh);
        auto row = m_walker_list.m_row;
        auto large_enough = [&](const size_t& isubject){return std::abs(row.m_weight[isubject]) > thresh*nadd;};
        for (row.restart(); row.in_range(); row.step()){
            if (row.is_cleared()) continue;
            if (std::all_of(m_subject_iparts.cbegin(), m_subject_iparts.cend(), large_enough)) add_wf_row(row);
        }
        build_hamiltonian(ham);
    }


    void build_from_highest_weighted(const WalkerTable& list, size_t ndet_tot) {
        /*
        defs::inds row_inds;

        FermionOnvConnection conn(ref);
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