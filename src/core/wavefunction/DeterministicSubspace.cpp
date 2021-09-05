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

defs::inds DeterministicSubspace::make_iparts() {
    if (m_wf.nreplica() == 1) return {m_iroot};
    auto ipart = m_wf.m_format.flatten({m_iroot, 0});
    return {ipart, ipart + 1};
}

void DeterministicSubspace::make_rdm_contribs(Rdms &rdms, const Mbf &ref, const std::forward_list<size_t> &icol_list) {
    auto &row_local = m_local.m_row;
    auto &row_global = m_global.m_row;

    auto icol_it = icol_list.cbegin();
    for (; icol_it != icol_list.cend(); ++icol_it) {
        row_global.jump(*icol_it);
        if (row_global.m_mbf == ref) continue;
        if (m_wf.nreplica() == 2) {
            rdms.make_contribs(row_local.m_mbf, row_global.m_mbf,
                               row_local.m_weight[0] * row_global.m_weight[1]);
            rdms.make_contribs(row_local.m_mbf, row_global.m_mbf,
                               row_local.m_weight[1] * row_global.m_weight[0]);
        } else {
            rdms.make_contribs(row_local.m_mbf, row_global.m_mbf, row_local.m_weight[0] * row_global.m_weight[0]);
        }
    }
}

DeterministicSubspace::DeterministicSubspace(
        const fciqmc_config::Semistochastic &opts, Wavefunction &wf, size_t iroot) :
        Wavefunction::PartSharedRowSet<DeterministicDataRow>(
                wf, "semistochastic", {wf}, DeterministicDataRow::load_fn),
        m_opts(opts), m_wf(wf), m_iroot(iroot), m_iparts(make_iparts()) {}

void DeterministicSubspace::add_(WalkerTableRow &row) {
    base_t::add_(row.index());
    row.m_deterministic.set(m_iroot);
}

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
    size_t n_hconn = 0ul;
    size_t n_rdm_conn = 0ul;
    m_ham_matrix.resize(m_local.m_hwm);
    m_rdm_network.resize(m_local.m_hwm);
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        // loop over local subspace (H rows)
        auto &row_global = m_global.m_row;
        for (row_global.restart(); row_global.in_range(); row_global.step()) {
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            conn_work.connect(row_local.m_mbf, row_global.m_mbf);
            const auto exsig = conn_work.exsig();
            if (!exsig || exsig > defs::nexsig) continue; // diagonal
            auto helem = ham.get_element(row_local.m_mbf, conn_work);
            if (!consts::float_is_zero(helem)) {
                m_ham_matrix.add(row_local.index(), row_global.index(), helem);
                ++n_hconn;
            } else if (bilinears.m_rdms.takes_contribs_from(conn_work.exsig())) {
                m_rdm_network.add(row_local.index(), row_global.index());
                ++n_rdm_conn;
            }
        }
    }
    log::info("Number of H-connected pairs of MBFs in the deterministic subspace: {}", n_hconn);
    if (bilinears.m_rdms)
        log::info("Number of H-unconnected, but RDM-contributing pairs of MBFs: {}", n_rdm_conn);
}

void DeterministicSubspace::make_rdm_contribs(Rdms &rdms, const field::Mbf &ref) {
    if (!rdms || !rdms.m_accum_epoch) return;
    auto &row_local = m_local.m_row;
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        if (row_local.m_mbf == ref) continue;
        make_rdm_contribs(rdms, ref, m_ham_matrix[row_local.index()].first);
        make_rdm_contribs(rdms, ref, m_rdm_network[row_local.index()]);
    }
}

void DeterministicSubspace::project(double tau) {
    auto &row_local = m_local.m_row;
    auto &row_global = m_global.m_row;
    auto &row_wf = m_wf.m_store.m_row;
    auto irow_wf_it = m_irows.cbegin();
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        row_wf.jump(*irow_wf_it);
        auto lists = m_ham_matrix[row_local.index()];
        auto icol_it = lists.first.cbegin();
        auto value_it = lists.second.cbegin();
        for (; icol_it != lists.first.cend(); (++icol_it, ++value_it)) {
            DEBUG_ASSERT_FALSE(value_it == lists.second.cend(), "values list incongruent with column indices list");
            row_global.jump(*icol_it);
            // one replica or two
            for (const auto &ipart: m_iparts)
                m_wf.change_weight(ipart, -*value_it * tau * row_global.m_weight[ipart]);
        }
        ++irow_wf_it;
    }
}

DeterministicSubspaces::DeterministicSubspaces(const fciqmc_config::Semistochastic &opts) :
        m_opts(opts), m_epoch("semistochastic") {
}

DeterministicSubspaces::operator bool() const {
    return m_opts.m_size && m_epoch;
}

void
DeterministicSubspaces::build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears, Wavefunction &wf,
                                                 size_t icycle) {
    m_detsubs.resize(wf.nroot());
    REQUIRE_FALSE_ALL(bool(*this), "epoch should not be started when building deterministic subspaces");
    for (size_t iroot = 0ul; iroot < wf.nroot(); ++iroot) {
        REQUIRE_TRUE_ALL(m_detsubs[iroot] == nullptr, "detsubs should not already be allocated");
        m_detsubs[iroot] = mem_utils::make_unique<DeterministicSubspace>(m_opts, wf, iroot);
        m_detsubs[iroot]->build_from_most_occupied(ham, bilinears);
    }
    m_epoch.update(icycle, true);
}

void DeterministicSubspaces::update() {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->update();
}

void DeterministicSubspaces::project(double tau) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->project(tau);
}

void DeterministicSubspaces::make_rdm_contribs(Rdms &rdms, const Mbf &ref) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->make_rdm_contribs(rdms, ref);
}
