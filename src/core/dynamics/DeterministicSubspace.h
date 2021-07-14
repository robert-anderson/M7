//
// Created by rja on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H

#include <src/core/observables/MevGroup.h>
#include <src/core/field/Onv.h>
#include "src/core/sparse/SparseMatrix.h"
#include "src/core/hamiltonian/FermionHamiltonian.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/Reducible.h"
#include "WalkerTable.h"
#include "Wavefunction.h"

/**
 * captures only the data from WalkerTableRow relevant to semistochastic propagation
 */
struct DeterministicDataRow : Row {
    fields::Onv<> m_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;

    fields::Onv<> &key_field();

    DeterministicDataRow(const Wavefunction& wf);

    static void load_fn(const WalkerTableRow& source, DeterministicDataRow& local);
};


class PureConnections {
    bool m_resized_by_add = false;
    std::vector<std::forward_list<size_t>> m_rows{};

public:
    void resize(const size_t nrow) {
        ASSERT(nrow >= m_rows.size());
        m_rows.resize(nrow);
    }

    void expand(const size_t delta_nrow) {
        resize(m_rows.size() + delta_nrow);
    }

    size_t nrow() const {
        return m_rows.size();
    }

    void add(const size_t &irow, const size_t &icol) {
        if (irow >= m_rows.size()) {
            if (!m_resized_by_add) {
                log::warn("Resizing SparseMatrix by adding a row (this entails reallocation which is inefficient)");
                log::warn("Call resize before filling if number of rows is known in advance");
                m_resized_by_add = true;
            }
            resize(irow + 1);
        }
        m_rows[irow].push_front(icol);
    }

    bool empty() { return m_rows.empty(); }

    const std::forward_list<size_t>& row(const size_t &irow) const {
        return m_rows[irow];
    }
};


struct DeterministicSubspace : Wavefunction::PartSharedRowSet<DeterministicDataRow>{
    typedef Wavefunction::PartSharedRowSet<DeterministicDataRow> base_t;
    const fciqmc_config::Semistochastic& m_opts;
    Wavefunction& m_wf;
    sparse::Matrix<defs::ham_t> m_sparse_ham;
    Epoch m_epoch;

    DeterministicSubspace(const fciqmc_config::Semistochastic& opts, Wavefunction& wf, size_t icycle);

    virtual ~DeterministicSubspace(){}

    void add_(WalkerTableRow& row) {
        base_t::add_(row.index());
        for (size_t ipart=0ul; ipart<m_wf.npart(); ++ipart)
            row.m_deterministic.set(ipart);
    }

    void build_from_most_occupied(const FermionHamiltonian &ham);

    void build_connections(const FermionHamiltonian &ham);

    void build_from_all_occupied(const FermionHamiltonian &ham);

    void build_from_occupied_connections(const FermionHamiltonian &ham, const fields::Onv<>& onv);

    void make_mev_contribs(MevGroup& mevs, const fields::Onv<>& ref);

    /**
      * for every deterministically-propagated row on this MPI rank, update its value.
      * @param tau
      *  timestep
      */
    void project(double tau);
};

#endif //M7_DETERMINISTICSUBSPACE_H