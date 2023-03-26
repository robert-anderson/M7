//
// Created by Robert J. Anderson on 16/06/2020.
//

#include "DeterministicSubspace.h"
#include "M7_lib/util/Pointer.h"


uintv_t deterministic::Subspace::make_iparts() {
    if (m_wf.nreplica() == 1) return {m_iroot};
    auto ipart = m_wf.m_format.flatten({m_iroot, 0});
    return {ipart, ipart + 1};
}

void deterministic::Subspace::make_rdm_contrib(Rdms &rdms, const shared_rows::Walker *hf, const sparse::Element& elem) {
    const auto& row = this->gathered().m_row;
    row.jump(elem);
    if (hf && (row.m_mbf == hf->mbf())) return;
    if (m_wf.nreplica() == 2) {
        rdms.make_contribs(m_local_row.m_mbf, row.m_mbf,
                           m_local_row.m_weight[0] * row.m_weight[1]);
        rdms.make_contribs(m_local_row.m_mbf, row.m_mbf,
                           m_local_row.m_weight[1] * row.m_weight[0]);
    } else {
        rdms.make_contribs(m_local_row.m_mbf, row.m_mbf, m_local_row.m_weight[0] * row.m_weight[0]);
    }
}

deterministic::Subspace::Subspace(
        const conf::Semistochastic &opts, wf::Vectors &wf, uint_t iroot) :
        shared_rows::Set<Walker>("semistochastic", wf.m_store),
        m_opts(opts), m_wf(wf),
        m_iroot(iroot), m_local_row(wf.m_store.m_row), m_iparts(make_iparts()){}

void deterministic::Subspace::add_(Walker &row) {
    base_t::add_(row.index());
    row.m_deterministic.set(m_iroot);
    logging::debug_("adding MBF {} to the deterministic subspace", row.m_mbf.to_string());
}

void deterministic::Subspace::select_highest_weighted() {
    auto row1 = m_wf.m_store.m_row;
    auto row2 = row1;
    wf::Vectors::weights_gxr_t gxr(row1.m_weight, row2.m_weight, true, true, m_iparts);
    gxr.find(m_opts.m_size);
    for (uint_t i = 0ul; i < gxr.m_ninclude.m_local; ++i) {
        row1.jump(gxr[i]);
        add_(row1);
    }
}

void deterministic::Subspace::select_l1_norm_fraction() {
    auto av_l1_norm = m_wf.m_stats.m_nwalker.total().sum_over(m_iparts) / m_iparts.size();
    const auto cutoff = av_l1_norm * m_opts.m_l1_fraction_cutoff.m_value;
    auto row = m_wf.m_store.m_row;
    for (row.restart(); row; ++row){
        av_l1_norm = row.m_weight.sum_over(m_iparts) / m_iparts.size();
        if (std::abs(av_l1_norm) >= cutoff) add_(row);
    }
}

void deterministic::Subspace::make_connections(const Hamiltonian &ham, const Rdms &rdms){
    const auto& gathered = this->gathered();
    full_update();
    if (!gathered.nrow_in_use()) {
        logging::info("No MBFs were selected to form a deterministic subspace for root {}", m_iroot);
        return;
    }
    logging::info("Forming a deterministic subspace with {} MBFs", gathered.nrow_in_use());
    suite::Conns conns_work(m_wf.m_sector.size());
    auto &conn_work = conns_work[m_local_row.m_mbf];
    uint_t nconn_ham = 0ul;
    uint_t nconn_rdm = 0ul;
    m_ham_matrix.resize(nrec_());
    m_rdm_network.resize(nrec_());
    uint_t iirow = ~0ul;
    for (auto irec: m_irecs) {
        ++iirow;
        m_local_row.jump(irec);
        // loop over local subspace (H rows)
        auto& all_row = gathered.m_row;
        for (all_row.restart(); all_row; ++all_row) {
            // loop over full subspace (H columns)
            // only add to sparse H if dets are connected
            if (mbf::exsig(m_local_row.m_mbf, all_row.m_mbf) == opsig::c_invalid) continue;
            conn_work.connect(m_local_row.m_mbf, all_row.m_mbf);
            const auto exsig = conn_work.exsig();
            if (exsig == opsig::c_zero) continue; // diagonal
            auto helem = ham.get_element(m_local_row.m_mbf, conn_work);
            if (ham::is_significant(helem)) {
                m_ham_matrix.add(iirow, {all_row.index(), helem});
                ++nconn_ham;
            } else if (rdms.takes_contribs_from(conn_work.exsig())) {
                m_rdm_network.add(iirow, all_row.index());
                ++nconn_rdm;
            }
        }
    }
    logging::info_("Size of local deterministic subspace: {}", nrec_());
    nconn_ham = mpi::all_sum(nconn_ham);
    nconn_rdm = mpi::all_sum(nconn_rdm);
    logging::info("Number of H-connected pairs of MBFs in the deterministic subspace: {}", nconn_ham);
    if (rdms) logging::info("Number of H-unconnected, but RDM-contributing pairs of MBFs: {}", nconn_rdm);
}

void deterministic::Subspace::make_connections(const Hamiltonian &ham, const SpecMoms &spec_moms) {
    if (!spec_moms) return;
    m_frm_hole_perturbed = ptr::smart::make_unique<FrmOpPerturbed>(
            m_wf.m_sector, true, spec_moms.m_selected_spinorbs, m_iroot);
    m_frm_hole_perturbed->setup(ham, gathered());
    m_frm_particle_perturbed = ptr::smart::make_unique<FrmOpPerturbed>(
            m_wf.m_sector, false, spec_moms.m_selected_spinorbs, m_iroot);
    m_frm_particle_perturbed->setup(ham, gathered());
}

void deterministic::Subspace::make_rdm_contribs(Rdms &rdms, const shared_rows::Walker *hf) {
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

void deterministic::Subspace::project(double tau) {
    // all required communication has already been done: the entire semistochastic vector is available on all ranks
    auto &row = gathered().m_row;
    uint_t iirec = ~0ul;
    // loop over row indices in the portion of the wavefunction stored on this rank
    for (auto irec: m_irecs) {
        ++iirec;
        auto& walker = m_wf.m_store.m_row;
        walker.jump(irec);
        v_t<wf_t> delta(m_wf.npart(), 0.0);
        for (const auto& elem: m_ham_matrix[iirec]){
            row.jump(elem.m_i);
            // one replica or two
            for (const auto &ipart: m_iparts) {
                // update the walker population locally due to coeffs from all connected deterministic MBFs
                const auto coeff = row.m_weight[ipart];
                delta[ipart] -= tau * elem.m_v * coeff;
            }
        }
        for (const auto &ipart: m_iparts) m_wf.change_weight(walker, ipart, delta[ipart]);
    }
}

deterministic::Subspaces::Subspaces(const conf::Semistochastic &opts) :
        m_opts(opts), m_epoch("semistochastic") {
}

deterministic::Subspaces::operator bool() const {
    return m_opts.m_enabled && m_epoch;
}

void deterministic::Subspaces::init(const Hamiltonian &ham, const Maes& maes, wf::Vectors &wf, uint_t icycle) {
    m_detsubs.resize(wf.nroot());
    REQUIRE_FALSE_ALL(bool(*this), "epoch should not be started when building deterministic subspaces");

    for (uint_t iroot = 0ul; iroot < wf.nroot(); ++iroot) {
        auto& detsub = m_detsubs[iroot];
        REQUIRE_TRUE_ALL(detsub == nullptr, "detsubs should not already be allocated");
        detsub = ptr::smart::make_unique<Subspace>(m_opts, wf, iroot);
        if (m_opts.m_l1_fraction_cutoff.m_value < 1.0) {
            logging::info("Selecting walkers with magnitude >= {:.5f}% of the current global population "
                          "for root {} deterministic subspace", 100.0*m_opts.m_l1_fraction_cutoff.m_value, iroot);
            detsub->select_l1_norm_fraction();
        } else {
            logging::info("Selecting upto {} largest-magnitude walkers for root {} deterministic subspace",
                          m_opts.m_size, iroot);
            detsub->select_highest_weighted();
        }
        detsub->make_connections(ham, maes.m_rdms);
        detsub->make_connections(ham, maes.m_spec_moms);
    }

    if (m_opts.m_save.m_enabled) {
        // the subspaces are not going to change, so might as well dump them to archive now
        hdf5::FileWriter fw(m_opts.m_save.m_path);
        for (auto& detsub: m_detsubs) detsub->save(fw);
    }

    m_epoch.update(icycle, true);
}

void deterministic::Subspaces::update() {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->update();
}

void deterministic::Subspaces::project(double tau) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->project(tau);
}

void deterministic::Subspaces::make_rdm_contribs(Rdms &rdms, const shared_rows::Walker *hf) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->make_rdm_contribs(rdms, hf);
}

void deterministic::Subspaces::make_spec_mom_contribs(SpecMoms& spec_moms) {
    if (!*this) return;
    for (auto &detsub: m_detsubs) detsub->make_spec_mom_contribs(spec_moms);
}

deterministic::FrmOpPerturbed::FrmOpPerturbed(sys::Sector sector, bool hole, uintv_t select_pert_inds, uint_t iroot) :
        m_frm_basis(sector.m_frm.m_basis), m_ms2_conserve(sector.m_frm.m_elecs.m_ms2.conserve()), m_hole(hole),
        m_pert_basis_table(logging::format("Perturbed ({}) deterministic subspace for root {}",
                                           hole ? "hole" : "particle", iroot), {pert_basis_row_t(sector, "mbf")}),
        m_select_pert_inds(std::move(select_pert_inds)){
        m_pert_basis_table.set_expansion_factor(2.0);
    /*
     * initialise map entries to empty
     */
    for (auto ispinorb: m_select_pert_inds) m_basis_maps.insert({ispinorb, {}});
}

bool deterministic::FrmOpPerturbed::is_selected_pertuber(uint_t ispinorb) const {
    auto it = std::find(m_select_pert_inds.cbegin(), m_select_pert_inds.cend(), ispinorb);
    return it != m_select_pert_inds.cend();
}

void deterministic::FrmOpPerturbed::setup_basis(const MappedTable<Walker>& walker_subspace, uint_t ispinorb) {
    // don't keep m_pert_basis_table rows for spin orbital indices which have not been selected as pertubers
    const auto is_selected = is_selected_pertuber(ispinorb);
    auto& basis_map_row = m_basis_maps[ispinorb];
    // need to modify MBF, so copy into a working object
    buffered::Mbf work_mbf(walker_subspace.m_row.m_mbf);
    conn::Mbf work_conn(work_mbf);
    // row in the walker deterministic subspace
    auto walker = walker_subspace.m_row;
    for (walker.restart(); walker; ++walker) {
        /*
         * ignore this N-electron state if it would be destroyed by application of the perturber:
         * if hole, ispinorb must be occupied, else, ispinorb must be vacant
         */
        if (mbf::get_spinorb(walker.m_mbf, ispinorb) != m_hole) continue;
        /*
         * if hole, put 0, else, put 1
         */
        set_conn(work_conn, ispinorb, m_hole);
        work_conn.apply(walker.m_mbf, work_mbf);
        const bool phase = work_conn.phase(walker.m_mbf);
        // index in the walker deterministic subspace
        const auto ici = walker.index();
        // see if this perturber has already been generated
        auto& lookup = m_pert_basis_table.lookup(work_mbf);
        const auto ipert = lookup ? lookup.index() : m_pert_basis_table.insert(work_mbf).index();
        if (is_selected) basis_map_row.push_back({ici, ipert, phase});
    }
    // m_pert_basis_table now contains the union of all perturbers
    // m_basis_map now contains correspondences between walker deterministic subspace and perturber basis
}

void deterministic::FrmOpPerturbed::setup_ham(const Hamiltonian& h) {
    // share out the perturbed-space Hamiltonian rows evenly
    const auto count_local = mpi::evenly_shared_count(full_basis_size());
    m_pert_basis_displ = mpi::evenly_shared_displ(full_basis_size());
    auto row = m_pert_basis_table.m_row;
    const auto& src = row.m_field;
    for (row.jump(m_pert_basis_displ); row.in_range(m_pert_basis_displ + count_local); ++row) {
        const auto helem_diag = h.get_element(src);
        const auto irow = row.index() - m_pert_basis_displ;
        DEBUG_ASSERT_TRUE(m_ham_pert[irow].empty(), "sparse Hamiltonian row should be empty");
        if (ham::is_significant(helem_diag)) m_ham_pert.insert(irow, {row.index(), helem_diag});
        auto& col = m_pert_basis_table.m_row;
        const auto& dst = col.m_field;
        for (col.restart(); col; ++col) {
            const auto icol = col.index();
            const auto helem = h.get_element(src, dst);
            if (ham::is_significant(helem)) m_ham_pert.insert(irow, {icol, helem});
        }
    }
    // m_perm_ham now contains the Hamiltonian projected into the perturbed space
}

void deterministic::FrmOpPerturbed::fill_full_vec(const MappedTable<Walker>& walker_subspace, uint_t ispinorb_right,
                                                  uint_t ipart_right) {
    auto& walker = walker_subspace.m_row;
    m_full_work_vec.assign(full_basis_size(), 0.0);
    for (auto& entry : m_basis_maps[ispinorb_right]) {
        walker.jump(entry.m_iwalker);
        m_full_work_vec[entry.m_iperturbed] += (entry.m_phase ? -1 : 1) * walker.m_weight[ipart_right];
    }
    // m_full_work_vec now contains the perturbed walker deterministic subspace with perturber ispinorb
}

wf_t deterministic::FrmOpPerturbed::contract(const MappedTable<Walker>& walker_subspace, uint_t ispinorb_left,
                                             uint_t ipart_left) {
    wf_t inner_product = 0.0;
    auto& walker = walker_subspace.m_row;
    for (auto& entry : m_basis_maps[ispinorb_left]) {
        const auto full_vec_elem = m_full_work_vec[entry.m_iperturbed];
        walker.jump(entry.m_iwalker);
        inner_product += (entry.m_phase ? -1 : 1) * full_vec_elem * walker.m_weight[ipart_left];
    }
    return inner_product;
}

void deterministic::FrmOpPerturbed::project_ham() {
    m_part_work_vec.assign(part_basis_size(), 0.0);
    for (uint_t irow = 0ul; irow < part_basis_size(); ++irow){
        const auto& row = m_ham_pert.get(irow);
        for (auto& entry: row) {
            m_part_work_vec[irow] += m_full_work_vec[entry.m_i] * entry.m_v;
        }
    }
    // bring all the parts together, overwriting the full vector
    mpi::all_gatherv(m_part_work_vec, m_full_work_vec);
}

void deterministic::FrmOpPerturbed::make_contribs(
        SpecMoms& spec_moms, const MappedTable<Walker>& walker_subspace, uint_t ipart, uint_t ipart_replica) {
    for (auto ispinorb_right : m_select_pert_inds) {
        auto* spec_mom = m_hole ? spec_moms.m_hole_spec_moms.data() : spec_moms.m_particle_spec_moms.data();
        // get spin so we can skip non-conserving left perturbers H is Ms2 conserving
        const auto spin = m_frm_basis.ispin(ispinorb_right);
        fill_full_vec(walker_subspace, ispinorb_right, ipart);
        for (uint_t ih=0ul; ih <= spec_moms.m_max_order; ++ih) {
            for (auto ispinorb_left : m_select_pert_inds) {
                // use hermiticity: skip
                // if (ispinorb_left < ispinorb_right) continue;
                // skip if spin-conservation would be violated
                if (m_ms2_conserve && m_frm_basis.ispin(ispinorb_left)!=spin) continue;
                auto inner_product = contract(walker_subspace, ispinorb_left, ipart_replica);
                spec_mom->make_contrib(ispinorb_left, ispinorb_right, inner_product);
            }
            project_ham();
            ++spec_mom;
        }
    }
}