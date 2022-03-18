//
// Created by Robert John Anderson on 2020-04-03.
//

#include <basis/Suites.h>
#include "Wavefunction.h"

Wavefunction::Wavefunction(const fciqmc_config::Document &opts, BasisData bd) :
        Communicator<WalkerTableRow, SpawnTableRow, false>(
                "wavefunction",
                opts.m_propagator.m_nw_target,
                size_t(opts.m_propagator.m_nw_target * opts.m_propagator.m_tau_init),
                opts.m_wavefunction.m_buffers,
                opts.m_wavefunction.m_load_balancing,
                {
                        {bd, opts.m_wavefunction.m_nroot,
                         opts.m_av_ests.any_bilinears() ? 2ul:1ul, need_av_weights(opts)},
                        MappedTableBase::nbucket_guess(
                                opts.m_propagator.m_nw_target / mpi::nrank(),
                                opts.m_wavefunction.m_hash_mapping.m_remap_ratio),
                        opts.m_wavefunction.m_hash_mapping.m_remap_nlookup,
                        opts.m_wavefunction.m_hash_mapping.m_remap_ratio
                },
                {{bd, need_send_parents(opts)}}),
        Archivable("wavefunction", opts.m_wavefunction.m_archivable),
        m_opts(opts),
        m_bd(bd),
        m_format(m_store.m_row.m_weight.m_format),
        m_ninitiator(m_format),
        m_delta_ninitiator(m_format),
        m_nwalker(m_format),
        m_delta_nwalker(m_format),
        m_l2_norm_square(m_format),
        m_delta_l2_norm_square(m_format),
        m_nspawned(m_format),
        m_nannihilated(m_format) {
    ASSERT(m_comm.recv().m_row.m_dst_mbf.belongs_to_row());
    m_summables.add_members(m_ninitiator, m_delta_ninitiator, m_nocc_mbf, m_delta_nocc_mbf,
                            m_nwalker, m_delta_nwalker, m_l2_norm_square, m_delta_l2_norm_square,
                            m_nspawned, m_nannihilated);
}

void Wavefunction::log_top_weighted(size_t ipart, size_t nrow) {
    weights_gxr_t gxr(m_store.m_row, m_store.m_row.m_weight, true, true, ipart);
    gxr.find(nrow);
    BufferedTable<WalkerTableRow> xr_gathered("global top weighted", {m_store.m_row});
    gxr.gatherv(xr_gathered);
    auto &row = xr_gathered.m_row;
    log::info("Top-weighted WF elements for part {}:", ipart);
    for (row.restart(); row.in_range(); row.step()) {
        log::info("{:<4} {}  {: .8e}  {}", row.index(), row.m_mbf, row.m_weight[ipart], row.m_initiator[ipart]);
    }
}

Wavefunction::~Wavefunction() {
    for (size_t ipart=0ul; ipart<npart(); ++ipart) log_top_weighted(ipart);
}

std::vector<std::string> Wavefunction::h5_field_names() {
    if (!defs::enable_bosons)
        return {"mbf", "weight"};
    else
        return {"mbf (fermion)", "mbf (boson)", "weight"};
}

void Wavefunction::h5_write(hdf5::GroupWriter &parent, std::string name) {
    m_store.save(parent, name, h5_field_names());
}

void Wavefunction::h5_read(hdf5::GroupReader &parent, const Hamiltonian &ham, const field::Mbf &ref,
                           std::string name) {
    m_store.clear();
    BufferedTable<WalkerTableRow> m_buffer("", {m_store.m_row});
    m_buffer.push_back();
    RowHdf5Reader<WalkerTableRow> row_reader(m_buffer.m_row, parent, name, h5_field_names());
    suite::Conns conn(m_bd);

    row_reader.restart();
    DEBUG_ASSERT_EQ(row_reader.m_weight.nelement(), m_format.m_nelement, "row reader has incompatible dimensionality");
    for (size_t iitem = 0ul; iitem < row_reader.m_nitem; ++iitem) {
        row_reader.read(iitem);
        conn[ref].connect(ref, row_reader.m_mbf);
        bool ref_conn = !consts::nearly_zero(ham.get_element(ref, conn[ref]));
        create_row(0ul, row_reader.m_mbf, ham.get_energy(row_reader.m_mbf), std::vector<bool>(npart(), ref_conn));
        set_weight(row_reader.m_weight);
    }
}

void Wavefunction::begin_cycle() {
    m_summables.zero_all_local();
    m_store.attempt_remap();
}

void Wavefunction::end_cycle() {
    m_summables.all_sum();
}

defs::wf_comp_t Wavefunction::square_norm(const size_t &ipart) const {
    defs::wf_comp_t res = 0.0;
    auto &row = m_store.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        const defs::wf_t &weight = row.m_weight[ipart];
        res += std::pow(std::abs(weight), 2.0);
    }
    return mpi::all_sum(res);
}

defs::wf_comp_t Wavefunction::l1_norm(const size_t &ipart) const {
    defs::wf_comp_t res = 0.0;
    auto &row = m_store.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        const defs::wf_t &weight = row.m_weight[ipart];
        res += std::abs(weight);
    }
    return mpi::all_sum(res);
}

void Wavefunction::grant_initiator_status(const size_t &ipart) {
    auto &row = m_store.m_row;
    if (row.m_initiator.get(ipart)) return;
    row.m_initiator.set(ipart);
    m_delta_ninitiator.m_local[ipart]++;
}

void Wavefunction::revoke_initiator_status(const size_t &ipart) {
    auto &row = m_store.m_row;
    if (!row.m_initiator.get(ipart)) return;
    row.m_initiator.clr(ipart);
    m_delta_ninitiator.m_local[ipart]--;
}

void Wavefunction::set_weight(const size_t &ipart, const defs::wf_t &new_weight) {
    auto &row = m_store.m_row;
    defs::wf_t &weight = row.m_weight[ipart];
    m_delta_nwalker.m_local[ipart] += std::abs(new_weight);
    m_delta_nwalker.m_local[ipart] -= std::abs(weight);
    m_delta_l2_norm_square.m_local[ipart] += std::pow(std::abs(new_weight), 2.0);
    m_delta_l2_norm_square.m_local[ipart] -= std::pow(std::abs(weight), 2.0);
    weight = new_weight;

    if (std::abs(new_weight) >= m_opts.m_propagator.m_nadd) grant_initiator_status(ipart);
    else revoke_initiator_status(ipart);
}

void Wavefunction::change_weight(const size_t &ipart, const defs::wf_t &delta) {
    set_weight(ipart, m_store.m_row.m_weight[ipart] + delta);
}

void Wavefunction::scale_weight(const size_t &ipart, const double &factor) {
    set_weight(ipart, factor * m_store.m_row.m_weight[ipart]);
}

void Wavefunction::zero_weight(const size_t &ipart) {
    set_weight(ipart, 0.0);
}

void Wavefunction::remove_row() {
    auto lookup = m_store[m_store.m_row.m_mbf];
    ASSERT(lookup);
    for (size_t ipart = 0ul; ipart < m_format.m_nelement; ++ipart) {
        zero_weight(ipart);
        // in the case that nadd==0.0, the set_weight method won't revoke:
        revoke_initiator_status(ipart);
        m_delta_nocc_mbf.m_local--;
    }
    m_store.erase(lookup);
}

size_t Wavefunction::create_row_(const size_t &icycle, const Mbf &mbf, const defs::ham_comp_t &hdiag,
                                 const std::vector<bool> &refconns) {
    DEBUG_ASSERT_EQ(refconns.size(), npart(), "should have as many reference rows as WF parts");
    DEBUG_ASSERT_TRUE(mpi::i_am(m_ra.get_rank(mbf)),
                      "this method should only be called on the rank responsible for storing the MBF");
    auto irow = m_store.insert(mbf);
    m_delta_nocc_mbf.m_local++;
    m_store.m_row.jump(irow);
    DEBUG_ASSERT_EQ(m_store.m_row.key_field(), mbf, "MBF was not properly copied into key field of WF row");
    m_store.m_row.m_hdiag = hdiag;
    for (size_t ipart=0ul; ipart < npart(); ++ipart)
        m_store.m_row.m_ref_conn.put(ipart, refconns[ipart]);
    /*
     * we need to be very careful here of off-by-one-like mistakes. the initial walker is "created" at the beginning
     * of MC cycle 0, and so the stats line output for cycle 0 will show that the number of walkers is the initial
     * occupation of the initial row. if a spawning event leads to the creation of another row, it is created on
     * iteration 1 even though it is added in the annihilating call of iteration 0. so, if this method is called in
     * the annihilating process of MC cycle i, it actually "becomes occupied" on cycle i+1.
     */
    if (storing_av_weights()) {
        m_store.m_row.m_icycle_occ = icycle+1;
        m_store.m_row.m_average_weight = 0;
    }
    return irow;
}

TableBase::Loc Wavefunction::create_row(const size_t &icycle, const Mbf &mbf, const defs::ham_comp_t &hdiag,
                                        const std::vector<bool> &refconns) {
    size_t irank = m_ra.get_rank(mbf);
    size_t irow;
    if (mpi::i_am(irank)) {
        irow = create_row_(icycle, mbf, hdiag, refconns);
    }
    mpi::bcast(irow, irank);
    return {irank, irow};
}

size_t Wavefunction::add_spawn(const field::Mbf &dst_mbf, const defs::wf_t &delta,
                               bool initiator, bool deterministic, size_t dst_ipart) {
    auto irank = m_ra.get_rank(dst_mbf);
    auto &dst_table = send(irank);

    auto &row = dst_table.m_row;
    row.push_back_jump();

    row.m_dst_mbf = dst_mbf;
    row.m_delta_weight = delta;
    row.m_src_initiator = initiator;
    row.m_src_deterministic = deterministic;
    row.m_ipart_dst = dst_ipart;
    return row.index();
}

size_t Wavefunction::add_spawn(const field::Mbf &dst_mbf, const defs::wf_t &delta, bool initiator, bool deterministic,
                               size_t dst_ipart, const field::Mbf &src_mbf, const defs::wf_t &src_weight) {
    auto irow = add_spawn(dst_mbf, delta, initiator, deterministic, dst_ipart);
    auto irank = m_ra.get_rank(dst_mbf);
    auto &row = send(irank).m_row;
    if (row.m_send_parents) {
        row.m_src_mbf = src_mbf;
        row.m_src_weight = src_weight;
    }
    DEBUG_ASSERT_NE(dst_mbf, src_mbf, "spawning diagonally");
    return irow;
}

void Wavefunction::load_fn(hdf5::GroupReader &parent) {

}

void Wavefunction::save_fn(hdf5::GroupWriter &parent) {

}