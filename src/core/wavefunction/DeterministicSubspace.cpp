//
// Created by rja on 16/06/2020.
//

#include "DeterministicSubspace.h"

field::Mbf &DeterministicDataRow::key_field() {
    return m_mbf;
}

DeterministicDataRow::DeterministicDataRow(const Wavefunction &wf) :
        m_mbf(this, wf.m_store.m_row.m_mbf.nsite(), "mbf"),
        m_weight(this, wf.m_store.m_row.m_weight.m_format, "weight") {}

void DeterministicDataRow::load_fn(const WalkerTableRow &source, DeterministicDataRow &local) {
    local.m_mbf = source.m_mbf;
    local.m_weight = source.m_weight;
}

DeterministicSubspace::DeterministicSubspace(
        const fciqmc_config::Semistochastic &opts, Wavefunction &wf, size_t iroot) :
        Wavefunction::PartSharedRowSet<DeterministicDataRow>(
                wf, "semistochastic", {wf}, DeterministicDataRow::load_fn),
        m_opts(opts), m_wf(wf), m_iroot(iroot), m_iparts(make_iparts()) {}

void DeterministicSubspace::build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears) {
    auto row = m_wf.m_store.m_row;
    Wavefunction::weights_gxr_t gxr(row, row.m_weight, true, true, m_iparts);
    gxr.find(m_opts.m_size);
    for (size_t i = 0ul; i < gxr.m_ninclude.m_local; ++i) {
        row.jump(gxr[i]);
        add_(row);
    }
    build_connections(ham, bilinears);
}

void DeterministicSubspace::build_connections(const Hamiltonian &ham, const Bilinears &bilinears) {
    update();
    log::info("Forming a deterministic subspace with {} ONVs", m_global.m_hwm);
    suite::Conns conns_work(m_wf.m_nsite);
    auto &row_local = m_local.m_row;
    auto &conn_work = conns_work[row_local.m_mbf];
    m_ham_matrix.resize(m_global.m_hwm);
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        // loop over local subspace (H rows)
        auto &row_global = m_global.m_row;
        for (row_global.restart(); row_global.in_range(); row_global.step()) {
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            conn_work.connect(row_local.m_mbf, row_global.m_mbf);
            if (!conn_work.exsig()) continue; // diagonal
            auto helem = ham.get_element(row_local.m_mbf, conn_work);
            if (!consts::float_is_zero(helem))
                m_ham_matrix.add(row_local.index(), row_global.index(), helem);
            else if (bilinears.m_rdms.takes_contribs_from(conn_work.exsig())) {
                m_rdm_network.add(row_local.index(), row_global.index());
            }
        }
    }
}

void DeterministicSubspace::make_rdm_contribs(Rdms &rdms, const field::Mbf &ref) {
    if (!rdms || !rdms.m_accum_epoch) return;
    auto &row_local = m_local.m_row;
    auto &row_global = m_global.m_row;
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        if (row_local.m_mbf == ref) continue;
        auto icol_list = m_ham_matrix[row_local.index()].first;
        auto icol_it = icol_list.cbegin();
        for (; icol_it != icol_list.cend(); ++icol_it) {
            row_global.jump(*icol_it);
            if (row_global.m_mbf == ref) continue;
            if (m_wf.nreplica() == 2) {
                rdms.make_contribs(row_local.m_mbf, row_global.m_mbf,
                                   0.5 * row_local.m_weight[0] * row_global.m_weight[1]);
                rdms.make_contribs(row_local.m_mbf, row_global.m_mbf,
                                   0.5 * row_local.m_weight[1] * row_global.m_weight[0]);
            } else {
                rdms.make_contribs(row_local.m_mbf, row_global.m_mbf, row_local.m_weight[0] * row_global.m_weight[0]);
            }
        }
    }
}

void DeterministicSubspace::project(double tau) {
    auto &row_local = m_local.m_row;
    auto &row_global = m_global.m_row;
    auto row_wf = m_wf.m_store.m_row;
    auto irow_wf_it = m_irows.cbegin();
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        row_wf.jump(*irow_wf_it);
        auto lists = m_ham_matrix[row_local.index()];
        auto icol_it = lists.first.cbegin();
        auto value_it = lists.second.cbegin();
        for (; icol_it != lists.first.cend(); (++icol_it, ++value_it)) {
            DEBUG_ASSERT_FALSE(value_it == lists.second.cend(), "values list incongruent with column indices list");
            row_global.jump(*icol_it);
            row_wf.m_weight.sub_scaled(*value_it * tau, row_global.m_weight);
        }
        ++irow_wf_it;
    }
}
