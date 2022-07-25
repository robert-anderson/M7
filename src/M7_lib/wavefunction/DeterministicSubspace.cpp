//
// Created by Robert J. Anderson on 16/06/2020.
//

#include "DeterministicSubspace.h"
#include "M7_lib/util/SmartPtr.h"

#if 0
field::Mbf &DeterministicDataRow::key_field() {
    return m_mbf;
}

DeterministicDataRow::DeterministicDataRow(const Wavefunction &wf) :
        m_mbf(this, wf.m_sector, "mbf"),
        m_weight(this, wf.m_store.m_row.m_weight.m_format, "weight") {}

void DeterministicDataRow::load_fn(const WalkerTableRow &source, DeterministicDataRow &local) {
    local.m_mbf = source.m_mbf;
    local.m_weight = source.m_weight;
}
#endif

uintv_t DeterministicSubspace::make_iparts() {
    if (m_wf.nreplica() == 1) return {m_iroot};
    auto ipart = m_wf.m_format.flatten({m_iroot, 0});
    return {ipart, ipart + 1};
}

void DeterministicSubspace::make_rdm_contrib(Rdms &rdms, const Mbf &ref, const sparse::Element& elem) {
    m_all.m_row.jump(elem);
    if (m_all.m_row.m_mbf == ref) return;
    if (m_wf.nreplica() == 2) {
        rdms.make_contribs(m_local_row.m_mbf, m_all.m_row.m_mbf,
                           m_local_row.m_weight[0] * m_all.m_row.m_weight[1]);
        rdms.make_contribs(m_local_row.m_mbf, m_all.m_row.m_mbf,
                           m_local_row.m_weight[1] * m_all.m_row.m_weight[0]);
    } else {
        rdms.make_contribs(m_local_row.m_mbf, m_all.m_row.m_mbf, m_local_row.m_weight[0] * m_all.m_row.m_weight[0]);
    }
}

DeterministicSubspace::DeterministicSubspace(
        const conf::Semistochastic &opts, Wavefunction &wf, uint_t iroot) :
        shared_rows::Set<WalkerTableRow>("semistochastic", wf),
        m_opts(opts), m_wf(wf), m_iroot(iroot), m_local_row(wf.m_store.m_row), m_iparts(make_iparts()){}

void DeterministicSubspace::add_(WalkerTableRow &row) {
    base_t::add_(row.index());
    row.m_deterministic.set(m_iroot);
}

void DeterministicSubspace::build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears) {
    auto row = m_wf.m_store.m_row;
    Wavefunction::weights_gxr_t gxr(row, row.m_weight, true, true, m_iparts);
    gxr.find(m_opts.m_size);
    for (uint_t i = 0ul; i < gxr.m_ninclude.m_local; ++i) {
        row.jump(gxr[i]);
        add_(row);
    }
    build_connections(ham, bilinears);
}

void DeterministicSubspace::build_connections(const Hamiltonian &ham, const Bilinears &bilinears) {
    full_update();
    logging::info("Forming a deterministic subspace with {} ONVs", m_all.m_hwm);
    suite::Conns conns_work(m_wf.m_sector.size());
    auto &conn_work = conns_work[m_local_row.m_mbf];
    uint_t n_hconn = 0ul;
    uint_t n_rdm_conn = 0ul;
    m_ham_matrix.resize(nrow_());
    m_rdm_network.resize(nrow_());
    uint_t iirow = ~0ul;
    for (auto irow: m_irows) {
        ++iirow;
        m_local_row.jump(irow);
        // loop over local subspace (H rows)
        auto& all_row = m_all.m_row;
        for (all_row.restart(); all_row.in_range(); all_row.step()) {
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            conn_work.connect(m_local_row.m_mbf, all_row.m_mbf);
            const auto exsig = conn_work.exsig();
            if (!exsig || exsig > exsig::c_ndistinct) continue; // diagonal
            auto helem = ham.get_element(m_local_row.m_mbf, conn_work);
            if (ham::is_significant(helem)) {
                m_ham_matrix.add(m_local_row.index(), {all_row.index(), helem});
                ++n_hconn;
            } else if (bilinears.m_rdms.takes_contribs_from(conn_work.exsig())) {
                m_rdm_network.add(m_local_row.index(), all_row.index());
                ++n_rdm_conn;
            }
        }
    }
    logging::info("Number of H-connected pairs of MBFs in the deterministic subspace: {}", n_hconn);
    if (bilinears.m_rdms)
        logging::info("Number of H-unconnected, but RDM-contributing pairs of MBFs: {}", n_rdm_conn);
}

void DeterministicSubspace::make_rdm_contribs(Rdms &rdms, const field::Mbf &ref) {
    if (!rdms || !rdms.m_accum_epoch) return;
    uint_t iirow = ~0ul;
    for (auto irow: m_irows) {
        ++iirow;
        m_local_row.jump(irow);
        if (m_local_row.m_mbf == ref) continue;
        /*
         * make contributions due to hamiltonian connections
         */
        for (auto& elem : m_ham_matrix[iirow]){
            make_rdm_contrib(rdms, ref, elem);
        }
        /*
         * make contributions due to RDM-only connections
         */
        for (auto& elem : m_rdm_network[iirow]){
            make_rdm_contrib(rdms, ref, elem);
        }
    }
}

void DeterministicSubspace::project(double tau) {
    auto &all_row = m_all.m_row;
    auto &row_wf = m_wf.m_store.m_row;
    uint_t iirow = ~0ul;
    for (auto irow: m_irows) {
        ++iirow;
        row_wf.jump(irow);
        const auto& elems = m_ham_matrix[iirow];
        v_t<wf_t> delta(m_wf.npart(), 0.0);
        for (auto& elem: elems){
            all_row.jump(elem.m_i);
            // one replica or two
            for (const auto &ipart: m_iparts) delta[ipart] = -elem.m_v * tau * all_row.m_weight[ipart];
        }
        for (const auto &ipart: m_iparts) m_wf.change_weight(ipart, delta[ipart]);
    }
}

DeterministicSubspaces::DeterministicSubspaces(const conf::Semistochastic &opts) :
        m_opts(opts), m_epoch("semistochastic") {
}

DeterministicSubspaces::operator bool() const {
    return m_opts.m_size && m_epoch;
}

void DeterministicSubspaces::build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears,
                                                 Wavefunction &wf, uint_t icycle) {
    m_detsubs.resize(wf.nroot());
    REQUIRE_FALSE_ALL(bool(*this), "epoch should not be started when building deterministic subspaces");
    for (uint_t iroot = 0ul; iroot < wf.nroot(); ++iroot) {
        REQUIRE_TRUE_ALL(m_detsubs[iroot] == nullptr, "detsubs should not already be allocated");
        m_detsubs[iroot] = smart_ptr::make_unique<DeterministicSubspace>(m_opts, wf, iroot);
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
