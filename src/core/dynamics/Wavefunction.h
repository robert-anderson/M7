//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/io/Options.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include "src/core/table/Communicator.h"
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/field/Fields.h"
#include "src/core/sort/QuickSorter.h"

struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow> {

    const Options &m_opts;
    const size_t m_nsite;

    NdFormat<defs::ndim_wf> m_format;
    // current "part" i.e. the flat element of the format
    size_t m_ipart = 0ul;

    ReductionSyndicate m_summables;

    ReductionMember<size_t, defs::ndim_wf> m_ninitiator;
    ReductionMember<int, defs::ndim_wf> m_delta_ninitiator;
    ReductionMember<size_t, defs::ndim_wf> m_nocc_onv;
    ReductionMember<int, defs::ndim_wf> m_delta_nocc_onv;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_l2_norm_square;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_l2_norm_square;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_nannihilated;

    MappedTable<UniqueOnvRow> m_unique_recvd_onvs;
    MappedTable<OnvRow> m_parent_recvd_onvs;

    Wavefunction(const Options &opts, size_t nsite) :
            Communicator<WalkerTableRow, SpawnTableRow, false>(
                    "walker",
                    opts.walker_buffer_expansion_factor,
                    opts.nload_balance_block_per_rank*mpi::nrank(),
                    opts.load_balance_period,
                    {
                        WalkerTableRow(nsite, opts.nroot, opts.nreplica),
                        MappedTableBase::nbucket_guess(opts.nwalker_target / mpi::nrank(), 3)
                    },
                    {SpawnTableRow(nsite, opts.rdm_rank>0)},
                    opts.acceptable_load_imbalance
            ),
            m_opts(opts),
            m_nsite(nsite),
            m_format({opts.nroot, opts.nreplica}),
            m_ninitiator(m_summables, m_format),
            m_delta_ninitiator(m_summables, m_format),
            m_nocc_onv(m_summables, m_format),
            m_delta_nocc_onv(m_summables, m_format),
            m_nwalker(m_summables, m_format),
            m_delta_nwalker(m_summables, m_format),
            m_l2_norm_square(m_summables, m_format),
            m_delta_l2_norm_square(m_summables, m_format),
            m_nannihilated(m_summables, m_format),
            m_unique_recvd_onvs({}, 100),
            m_parent_recvd_onvs({nsite}, 100){
        m_store.resize((m_opts.walker_buffer_size_factor_initial*m_opts.nwalker_target)/mpi::nrank());
        m_comm.resize((m_opts.spawn_buffer_size_factor_initial*m_opts.nwalker_target)/mpi::nrank());
        ASSERT(m_comm.recv().m_row.m_dst_onv.is_added_to_row());
    }

    std::vector<std::string> h5_field_names() {
        if (!defs::enable_bosons)
            return {"onv", "weight"};
        else
            return {"onv (fermion)", "onv (boson)", "weight"};
    }

    void h5_write(hdf5::GroupWriter& parent, std::string name="wavefunction") {
        m_store.write(parent, name, h5_field_names());
    }

    void h5_read(hdf5::GroupReader& parent, const Hamiltonian<>& ham, const fields::Onv<>& ref, std::string name="wavefunction") {
        m_store.clear();
        BufferedTable<WalkerTableRow> m_buffer("", {{m_nsite, m_opts.nroot, m_opts.nreplica}});
        m_buffer.push_back();
        RowHdf5Reader<WalkerTableRow> row_reader(m_buffer.m_row, parent, name, h5_field_names());
        conn::Antisym<> conn(m_nsite);

        row_reader.restart();
        for (size_t iitem = 0ul; iitem<row_reader.m_nitem; ++iitem){
            row_reader.read(iitem);
            conn.connect(ref, row_reader.m_onv);
            bool ref_conn = conn.connected();
            conn.connect(row_reader.m_onv, row_reader.m_onv);
            create_walker_(row_reader.m_onv, row_reader.m_weight[0], ham.get_element(conn), ref_conn);
        }
    }


    void begin_cycle() {
        m_summables.zero();
        m_store.remap_if_required();
    }

    void end_cycle() {
        m_summables.all_sum();
    }

//    void update(size_t icycle, double work_time) {
//        m_ra.update(icycle, work_time, m_walkers, m_walkers.m_key_field);
//    }


    defs::wf_comp_t square_norm() const {
        defs::wf_comp_t res = 0.0;
        auto& row = m_store.m_row;
        for (row.restart(); row.in_range(); row.step()) {
            const defs::wf_t& weight = row.m_weight[m_ipart];
            res += std::pow(weight, 2.0);
        }
        return mpi::all_sum(res);
    }

    defs::wf_comp_t l1_norm() const {
        defs::wf_comp_t res = 0.0;
        auto& row = m_store.m_row;
        for (row.restart(); row.in_range(); row.step()) {
            const defs::wf_t& weight = row.m_weight[m_ipart];
            res += std::abs(weight);
        }
        return mpi::all_sum(res);
    }

    void grant_initiator_status() {
        auto& row = m_store.m_row;
        //MPI_ASSERT(!row.m_initiator.get(m_ipart), "row is already initiator");
        if (row.m_initiator.get(m_ipart)) return;
        row.m_initiator.set(m_ipart);
        m_delta_ninitiator(0, 0)++;
    }

    void revoke_initiator_status() {
        auto& row = m_store.m_row;
        //MPI_ASSERT(!row.m_initiator.get(m_ipart), "row is not initiator");
        if (!row.m_initiator.get(m_ipart)) return;
        row.m_initiator.clr(m_ipart);
        m_delta_ninitiator(0, 0)--;
    }

    void set_weight(const defs::wf_t &new_weight) {
        auto& row = m_store.m_row;
        defs::wf_t& weight = row.m_weight[m_ipart];
        m_delta_nwalker(0, 0) += std::abs(new_weight);
        m_delta_nwalker(0, 0) -= std::abs(weight);
        m_delta_l2_norm_square(0, 0) += std::pow(std::abs(new_weight), 2.0);
        m_delta_l2_norm_square(0, 0) -= std::pow(std::abs(weight), 2.0);
        weight = new_weight;

        if (std::abs(new_weight) >= m_opts.nadd_initiator) grant_initiator_status();
        else revoke_initiator_status();
    }

    void change_weight(const defs::wf_t &delta) {
        set_weight(m_store.m_row.m_weight[m_ipart] + delta);
    }

    void scale_weight(const double &factor) {
        set_weight(factor*m_store.m_row.m_weight[m_ipart]);
    }

    void zero_weight() {
        set_weight(0.0);
    }

    void remove_walker() {
        if (m_ra.row_mapped_by_dependent(m_store.m_row.m_i)) return;
        auto lookup = m_store[m_store.m_row.m_onv];
        ASSERT(lookup);
        zero_weight();
        // in the case that nadd==0.0, the set_weight method won't revoke:
        revoke_initiator_status();
        m_store.erase(lookup);
        m_delta_nocc_onv(0, 0)--;
    }

    size_t create_walker_(const fields::Onv<> &onv, const defs::ham_t weight,
                          const defs::ham_comp_t &hdiag, bool refconn) {
        ASSERT(mpi::i_am(m_ra.get_rank(onv)));
        if (m_store.is_full()) m_store.expand(1);
        auto irow = m_store.insert(onv);
        m_delta_nocc_onv(0, 0)++;
        m_store.m_row.jump(irow);
        ASSERT(m_store.m_row.m_onv == onv)
        set_weight(weight);
        m_store.m_row.m_hdiag = hdiag;
        m_store.m_row.m_reference_connection.put(m_ipart, refconn);
        m_store.m_row.m_deterministic.clr(m_ipart);
        return irow;
    }


    TableBase::Loc create_walker(const fields::Onv<> &onv, const defs::ham_t weight,
                                 const defs::ham_comp_t &hdiag, bool refconn) {
        size_t irank = m_ra.get_rank(onv);
        size_t irow;
        if (mpi::i_am(irank)) irow = create_walker_(onv, weight, hdiag, refconn);
        mpi::bcast(irow, irank);
        return {irank, irow};
    }

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart) {
        auto irank = m_ra.get_rank(dst_onv);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "SENDING SPAWNED WALKER" << std::endl;
        std::cout << consts::verb << "generated determinant:   " << dst_onv.to_string() << std::endl;
        std::cout << consts::verb << "destination rank:        " << irank << std::endl;
        std::cout << consts::verb << "spawned weight:          " << delta << std::endl;
        std::cout << consts::verb << "parent is initiator:     " << initiator << std::endl;
        std::cout << consts::verb << "parent is deterministic: " << deterministic << std::endl;
#endif
        auto &dst_table = send(irank);

        auto irow = dst_table.push_back();
        auto& row = dst_table.m_row;
        row.jump(irow);

        row.m_dst_onv = dst_onv;
        row.m_delta_weight = delta;
        row.m_src_initiator = initiator;
        row.m_src_deterministic = deterministic;
        row.m_dst_ipart = dst_ipart;
        return irow;
    }

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart,
                     const fields::Onv<> &src_onv, const defs::wf_t &src_weight) {
        auto irow = add_spawn(dst_onv, delta, initiator, deterministic, dst_ipart);
        auto irank = m_ra.get_rank(dst_onv);
        auto& row = send(irank).m_row;
        row.jump(irow);
        if (row.m_send_parents){
            row.m_src_onv = src_onv;
            row.m_src_weight = src_weight;
        }
        return irow;
    }

    void consolidate_spawned(){
//        auto row1 = recv().m_row;
//        auto row2 = recv().m_row;
//        auto comp_fn = [&](const size_t &irow1, const size_t &irow2){
//            row1.jump(irow1);
//            row2.jump(irow2);
//            return row1.m_dst_onv <= row2.m_dst_onv;
//        };
//        /*
//         * sorting in ascending lexical order
//         */
//        QuickSorter qs(comp_fn);
//        qs.sort(recv());
//
//        auto& parent_row = m_parent_recvd_onvs.m_row;
//        const auto &row = recv().m_row;
//        for (row.restart(); row.in_range(); row.step()){
//
//        }
    }
};

#endif //M7_WAVEFUNCTION_H
