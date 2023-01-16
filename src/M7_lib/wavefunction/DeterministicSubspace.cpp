//
// Created by Robert J. Anderson on 16/06/2020.
//

#include "DeterministicSubspace.h"
#include "M7_lib/util/Pointer.h"

#if 0
field::Mbf &DeterministicDataRow::key_field() {
    return m_mbf;
}

DeterministicDataRow::DeterministicDataRow(const Wavefunction &wf) :
        m_mbf(this, wf.m_sector, "mbf"),
        m_weight(this, wf.m_store.m_row.m_weight.m_format, "weight") {}

void DeterministicDataRow::load_fn(const Walker &source, DeterministicDataRow &local) {
    local.m_mbf = source.m_mbf;
    local.m_weight = source.m_weight;
}
#endif

uintv_t DeterministicSubspace::make_iparts() {
    if (m_wf.nreplica() == 1) return {m_iroot};
    auto ipart = m_wf.m_format.flatten({m_iroot, 0});
    return {ipart, ipart + 1};
}

void DeterministicSubspace::make_rdm_contrib(Rdms &rdms, const shared_rows::Walker *hf, const sparse::Element& elem) {
    m_all.m_row.jump(elem);
    if (hf && (m_all.m_row.m_mbf == hf->mbf())) return;
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
        shared_rows::Set<Walker>("semistochastic", wf.m_store),
        m_opts(opts), m_wf(wf), m_iroot(iroot), m_local_row(wf.m_store.m_row), m_iparts(make_iparts()){}

void DeterministicSubspace::add_(Walker &row) {
    base_t::add_(row.index());
    row.m_deterministic.set(m_iroot);
    logging::debug_("adding MBF {} to the deterministic subspace", row.m_mbf.to_string());
}

void DeterministicSubspace::select_highest_weighted() {
    auto row1 = m_wf.m_store.m_row;
    auto row2 = row1;
    Wavefunction::weights_gxr_t gxr(row1.m_weight, row2.m_weight, true, true, m_iparts);
    gxr.find(m_opts.m_size);
    for (uint_t i = 0ul; i < gxr.m_ninclude.m_local; ++i) {
        row1.jump(gxr[i]);
        add_(row1);
    }
}

void DeterministicSubspace::select_l1_norm_fraction() {
    auto av_l1_norm = m_wf.m_nwalker.m_reduced.sum_over(m_iparts) / m_iparts.size();
    const auto cutoff = av_l1_norm * m_opts.m_l1_fraction_cutoff.m_value;
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row; ++row){
        av_l1_norm = row.m_weight.sum_over(m_iparts) / m_iparts.size();
        if (av_l1_norm >= cutoff) add_(row);
    }
}

void DeterministicSubspace::make_connections(const Hamiltonian &ham, const Bilinears &bilinears) {
    full_update();
    logging::info("Forming a deterministic subspace with {} MBFs", m_all.nrow_in_use());
    suite::Conns conns_work(m_wf.m_sector.size());
    auto &conn_work = conns_work[m_local_row.m_mbf];
    uint_t n_hconn = 0ul;
    uint_t n_rdm_conn = 0ul;
    m_ham_matrix.resize(nrec_());
    m_rdm_network.resize(nrec_());
    uint_t iirow = ~0ul;
    for (auto irec: m_irecs) {
        ++iirow;
        m_local_row.jump(irec);
        // loop over local subspace (H rows)
        auto& all_row = m_all.m_row;
        for (all_row.restart(); all_row; ++all_row) {
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            conn_work.connect(m_local_row.m_mbf, all_row.m_mbf);
            const auto exsig = conn_work.exsig();
            if ((exsig == opsig::c_zero) || (exsig == opsig::c_invalid)) continue; // diagonal or out of bounds
            auto helem = ham.get_element(m_local_row.m_mbf, conn_work);
            if (ham::is_significant(helem)) {
                m_ham_matrix.add(iirow, {all_row.index(), helem});
                ++n_hconn;
            } else if (bilinears.m_rdms.takes_contribs_from(conn_work.exsig())) {
                m_rdm_network.add(iirow, all_row.index());
                ++n_rdm_conn;
            }
        }
    }
    logging::info("Number of H-connected pairs of MBFs in the deterministic subspace: {}", n_hconn);
    if (bilinears.m_rdms)
        logging::info("Number of H-unconnected, but RDM-contributing pairs of MBFs: {}", n_rdm_conn);
}

void DeterministicSubspace::make_rdm_contribs(Rdms &rdms, const shared_rows::Walker *hf) {
    if (!rdms || !rdms.m_accum_epoch) return;
    uint_t iirec = ~0ul;
    for (auto irec: m_irecs) {
        ++iirec;
        m_local_row.jump(irec);
        if (hf && (m_local_row.m_mbf == hf->mbf())) continue;
        /*
         * make contributions due to hamiltonian connections
         */
        for (auto& elem : m_ham_matrix[iirec]){
            make_rdm_contrib(rdms, hf, elem);
        }
        /*
         * make contributions due to RDM-only connections
         */
        for (auto& elem : m_rdm_network[iirec]){
            make_rdm_contrib(rdms, hf, elem);
        }
    }
}

void DeterministicSubspace::project(double tau) {
    auto &all_row = m_all.m_row;
    uint_t iirec = ~0ul;
    // loop over row indices in the portion of the wavefunction stored on this rank
    for (auto irec: m_irecs) {
        ++iirec;
        auto& walker = m_wf.m_store.m_row;
        walker.jump(irec);
        v_t<wf_t> delta(m_wf.npart(), 0.0);
        for (const auto& elem: m_ham_matrix[iirec]){
            all_row.jump(elem.m_i);
            // one replica or two
            for (const auto &ipart: m_iparts) {
                // update the walker population locally due to coeffs from all connected deterministic MBFs
                const auto coeff = all_row.m_weight[ipart];
                delta[ipart] -= tau * elem.m_v * coeff;
            }
        }
        for (const auto &ipart: m_iparts) m_wf.change_weight(walker, ipart, delta[ipart]);
    }
}

DeterministicSubspaces::DeterministicSubspaces(const conf::Semistochastic &opts) :
        m_opts(opts), m_epoch("semistochastic") {
}

DeterministicSubspaces::operator bool() const {
    return m_opts.m_size && m_epoch;
}

void DeterministicSubspaces::init(const Hamiltonian &ham, const Bilinears &bilinears,
                                  Wavefunction &wf, uint_t icycle) {
    m_detsubs.resize(wf.nroot());
    REQUIRE_FALSE_ALL(bool(*this), "epoch should not be started when building deterministic subspaces");


    for (uint_t iroot = 0ul; iroot < wf.nroot(); ++iroot) {
        auto& detsub = m_detsubs[iroot];
        REQUIRE_TRUE_ALL(detsub == nullptr, "detsubs should not already be allocated");
        detsub = ptr::smart::make_unique<DeterministicSubspace>(m_opts, wf, iroot);
        if (m_opts.m_l1_fraction_cutoff.m_value < 1.0) {
            logging::info("Selecting walkers with magnitude >= {:.2f}% of the current global population "
                          "for root {} deterministic subspace", m_opts.m_l1_fraction_cutoff, iroot);
            detsub->select_l1_norm_fraction();
        } else {
            logging::info("Selecting upto {} largest-magnitude walkers for root {} deterministic subspace",
                          m_opts.m_size, iroot);
            detsub->select_highest_weighted();
        }
        detsub->make_connections(ham, bilinears);
    }

    if (m_opts.m_save.m_enabled) {
        // the subspaces are not going to change, so might as well dump them to archive now
        hdf5::FileWriter fw(m_opts.m_save.m_path);
        for (auto& detsub: m_detsubs) detsub->save(fw);
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

void DeterministicSubspaces::make_rdm_contribs(Rdms &rdms, const shared_rows::Walker *hf) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->make_rdm_contribs(rdms, hf);
}
