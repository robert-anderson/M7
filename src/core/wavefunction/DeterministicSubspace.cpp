//
// Created by rja on 16/06/2020.
//

#include <src/core/basis/Suites.h>
#include "DeterministicSubspace.h"

fields::Onv<> &DeterministicDataRow::key_field() {
    return m_onv;
}

DeterministicDataRow::DeterministicDataRow(const Wavefunction &wf) :
        m_onv(this, onv::nsite(wf.m_store.m_row.m_onv), "onv"),
        m_weight(this, wf.m_store.m_row.m_weight.m_format, "weight"){}

void DeterministicDataRow::load_fn(const WalkerTableRow &source, DeterministicDataRow &local) {
    local.m_onv = source.m_onv;
    local.m_weight = source.m_weight;
}

DeterministicSubspace::DeterministicSubspace(const fciqmc_config::Semistochastic& opts, Wavefunction &wf, size_t icycle) :
        Wavefunction::PartSharedRowSet<DeterministicDataRow>(wf, "semistochastic", {wf}, DeterministicDataRow::load_fn),
        m_opts(opts), m_wf(wf), m_epoch("semistochastic"){
    m_epoch.update(icycle, true);
}

void DeterministicSubspace::build_from_most_occupied(const FermionHamiltonian &ham) {
    auto row = m_wf.m_store.m_row;
    Wavefunction::weights_gxr_t gxr(row, row.m_weight, true, true, 0);
    gxr.find(m_opts.m_size);
    for (size_t i =0ul; i<gxr.m_ninclude; ++i) {
        row.jump(gxr[i]);
        add_(row);
    }
    build_connections(ham);
}

void DeterministicSubspace::build_connections(const FermionHamiltonian &ham) {
    update();
    log::debug("Forming a deterministic subspace with {} ONVs", m_global.m_hwm);
    suite::Conns conns_work(m_wf.m_nsite);
    auto& row_local = m_local.m_row;
    auto& conn_work = conns_work[row_local.m_onv];
    m_sparse_ham.resize(m_global.m_hwm);
    for (row_local.restart(); row_local.in_range(); row_local.step()){
        // loop over local subspace (H rows)
        auto& row_global = m_global.m_row;
        for (row_global.restart(); row_global.in_range(); row_global.step()){
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            conn_work.connect(row_local.m_onv, row_global.m_onv);
            auto helem = ham.get_element(row_local.m_onv, conn_work);
            if (conn_work.size() > 0 && conn_work.size() < 5)
                m_sparse_ham.add(row_local.index(), row_global.index(), helem);
        }
    }
}

void DeterministicSubspace::build_from_all_occupied(const FermionHamiltonian &ham) {
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (!row.is_cleared()) add_(row);
    }
    build_connections(ham);
}

void DeterministicSubspace::build_from_occupied_connections(const FermionHamiltonian &ham, const fields::Onv<> &onv) {
    suite::Conns conns_work(m_wf.m_nsite);
    auto row = m_wf.m_store.m_row;
    auto& conn_work = conns_work[row.m_onv];
    for (row.restart(); row.in_range(); row.step()){
        conn_work.connect(onv, row.m_onv);
        if (row.is_cleared() || conn_work.size()>4) continue;
        add_(row);
        for (size_t ipart=0ul; ipart<m_wf.npart(); ++ipart)
            row.m_deterministic.set(ipart);
    }
    build_connections(ham);
}

void DeterministicSubspace::make_mev_contribs(MevGroup &mevs, const fields::Onv<> &ref) {
    if (!mevs.m_accum_epoch) return;
    auto& row_local = m_local.m_row;
    auto& row_global = m_global.m_row;
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        if (row_local.m_onv==ref) continue;
        for (const auto& entry: m_sparse_ham.row(row_local.index())) {
            row_global.jump(entry.icol);
            if (row_global.m_onv==ref) continue;
            mevs.m_fermion_rdm->make_contribs(row_local.m_onv, 0.5*row_local.m_weight[0], row_global.m_onv, row_global.m_weight[1]);
            mevs.m_fermion_rdm->make_contribs(row_local.m_onv, 0.5*row_local.m_weight[1], row_global.m_onv, row_global.m_weight[0]);
        }
    }
}

void DeterministicSubspace::project(double tau) {
    auto& row_local = m_local.m_row;
    auto& row_global = m_global.m_row;
    auto row_wf = m_wf.m_store.m_row;
    auto irow_wf_it = m_irows.cbegin();
    for (row_local.restart(); row_local.in_range(); row_local.step()) {
        row_wf.jump(*irow_wf_it);
        for (const auto& entry: m_sparse_ham.row(row_local.index())) {
            row_global.jump(entry.icol);
            row_wf.m_weight.sub_scaled(tau * entry.element, row_global.m_weight);
        }
        ++irow_wf_it;
    }
}
