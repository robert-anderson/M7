//
// Created by Robert John Anderson on 2020-04-03.
//

#include "Wavefunction.h"

Wavefunction::Wavefunction(const fciqmc_config::Document &opts, size_t nsite):
        Communicator<WalkerTableRow, SpawnTableRow, false>(
                "wavefunction",
                opts.m_propagator.m_nw_target,
                size_t(opts.m_propagator.m_nw_target*opts.m_propagator.m_tau_init),
                opts.m_wavefunction.m_buffers,
                opts.m_wavefunction.m_load_balancing,
                {
                        {
                                nsite, opts.m_wavefunction.m_nroot,
                                opts.m_wavefunction.m_replicate ? 2ul : 1ul,
                                need_av_weights(opts)
                        },
                        MappedTableBase::nbucket_guess(
                                opts.m_propagator.m_nw_target / mpi::nrank(),
                                opts.m_wavefunction.m_hash_mapping.m_remap_ratio),
                        opts.m_wavefunction.m_hash_mapping.m_remap_nlookup,
                        opts.m_wavefunction.m_hash_mapping.m_remap_ratio
                },
                {{nsite, need_send_parents(opts)}}),
        Archivable("wavefunction", opts.m_wavefunction.m_archivable),
        m_opts(opts),
        m_nsite(nsite),
        m_format(m_store.m_row.m_weight.m_format),
        m_ninitiator(m_format),
        m_delta_ninitiator(m_format),
        m_nwalker(m_format),
        m_delta_nwalker(m_format),
        m_l2_norm_square(m_format),
        m_delta_l2_norm_square(m_format),
        m_nspawned(m_format),
        m_nannihilated(m_format){
    m_store.expand((m_opts.m_wavefunction.m_buffers.m_store_fac_init * m_opts.m_propagator.m_nw_target) / mpi::nrank());
    m_comm.expand((m_opts.m_wavefunction.m_buffers.m_comm_fac_init * m_opts.m_propagator.m_nw_target) / mpi::nrank());
    ASSERT(m_comm.recv().m_row.m_dst_onv.belongs_to_row());
    m_summables.add_members(m_ninitiator, m_delta_ninitiator, m_nocc_onv, m_delta_nocc_onv,
                            m_nwalker, m_delta_nwalker, m_l2_norm_square, m_delta_l2_norm_square,
                            m_nspawned, m_nannihilated);
}

Wavefunction::~Wavefunction() {
    weights_gxr_t gxr(m_store.m_row, m_store.m_row.m_weight, true, true, 0);
    gxr.find(20);
    BufferedTable<WalkerTableRow> xr_gathered("global top weighted", {m_store.m_row});
    gxr.gatherv(xr_gathered);
    auto& row = xr_gathered.m_row;
    log::info("Top-weighted WF elements for part 0:");
    for (row.restart(); row.in_range(); row.step()){
        log::info("{}  {}  {}", row.m_onv.to_string(), row.m_weight[0], row.m_initiator[0]);
    }
}

std::vector<std::string> Wavefunction::h5_field_names() {
    if (!defs::enable_bosons)
        return {"onv", "weight"};
    else
        return {"onv (fermion)", "onv (boson)", "weight"};
}

void Wavefunction::h5_write(hdf5::GroupWriter &parent, std::string name) {
    m_store.save(parent, name, h5_field_names());
}

void Wavefunction::h5_read(hdf5::GroupReader &parent, const Hamiltonian<> &ham, const fields::Onv<> &ref,
                           std::string name) {
    m_store.clear();
    BufferedTable<WalkerTableRow> m_buffer("", {m_store.m_row});
    m_buffer.push_back();
    RowHdf5Reader<WalkerTableRow> row_reader(m_buffer.m_row, parent, name, h5_field_names());
    conn::Antisym<> conn(m_nsite);

    row_reader.restart();
    for (size_t iitem = 0ul; iitem < row_reader.m_nitem; ++iitem) {
        row_reader.read(iitem);
        conn.connect(ref, row_reader.m_onv);
        bool ref_conn = !consts::float_is_zero(ham.get_element(conn));
        ASSERT(row_reader.m_weight.nelement()==m_format.m_nelement);
        create_row(0ul, row_reader.m_onv, ham.get_energy(row_reader.m_onv), std::vector<bool>(npart(), ref_conn));
        set_weight(row_reader.m_weight);
    }
}

void Wavefunction::begin_cycle() {
    m_summables.zero_all_local();
    m_store.attempt_remap();
}

void Wavefunction::end_cycle() {
    m_summables.all_sum();
    m_store.attempt_remap();
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
    auto lookup = m_store[m_store.m_row.m_onv];
    ASSERT(lookup);
    for (size_t ipart = 0ul; ipart<m_format.m_nelement; ++ipart) {
        zero_weight(ipart);
        // in the case that nadd==0.0, the set_weight method won't revoke:
        revoke_initiator_status(ipart);
        m_delta_nocc_onv.m_local--;
    }
    m_store.erase(lookup);
}


size_t Wavefunction::add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                               bool initiator, bool deterministic, size_t dst_ipart) {
    auto irank = m_ra.get_rank(dst_onv);
    auto &dst_table = send(irank);

    auto &row = dst_table.m_row;
    row.push_back_jump();

    row.m_dst_onv = dst_onv;
    row.m_delta_weight = delta;
    row.m_src_initiator = initiator;
    row.m_src_deterministic = deterministic;
    row.m_dst_ipart = dst_ipart;
    return row.index();
}

size_t Wavefunction::add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta, bool initiator, bool deterministic,
                        size_t dst_ipart, const fields::Onv<> &src_onv, const defs::wf_t &src_weight) {
    auto irow = add_spawn(dst_onv, delta, initiator, deterministic, dst_ipart);
    auto irank = m_ra.get_rank(dst_onv);
    auto &row = send(irank).m_row;
    if (row.m_send_parents) {
        row.m_src_onv = src_onv;
        row.m_src_weight = src_weight;
    }
    return irow;
}

void Wavefunction::load_fn(hdf5::GroupReader &parent) {

}

void Wavefunction::save_fn(hdf5::GroupWriter &parent) {

}
