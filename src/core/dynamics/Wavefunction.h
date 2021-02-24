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

struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow> {

    const Options &m_opts;

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

    Wavefunction(const Options &opts, size_t nsite) :
            Communicator<WalkerTableRow, SpawnTableRow>(
                    "walker",
                    opts.walker_buffer_expansion_factor,
                    opts.nload_balance_block_per_rank*mpi::nrank(),
                    opts.load_balance_period,
                    WalkerTableRow(nsite, opts.nroot, opts.nreplica), SpawnTableRow(nsite),
                    MappedTableBase::nbucket_guess(opts.nwalker_target / mpi::nrank(), 3),
                    opts.acceptable_load_imbalance
            ),
            m_opts(opts),
            m_format({opts.nroot, opts.nreplica}),
            m_ninitiator(m_summables, m_format),
            m_delta_ninitiator(m_summables, m_format),
            m_nocc_onv(m_summables, m_format),
            m_delta_nocc_onv(m_summables, m_format),
            m_nwalker(m_summables, m_format),
            m_delta_nwalker(m_summables, m_format),
            m_l2_norm_square(m_summables, m_format),
            m_delta_l2_norm_square(m_summables, m_format),
            m_nannihilated(m_summables, m_format) {
        m_store.resize((m_opts.walker_buffer_size_factor_initial*m_opts.nwalker_target)/mpi::nrank());
        m_comm.resize((m_opts.spawn_buffer_size_factor_initial*m_opts.nwalker_target)/mpi::nrank());
        ASSERT(m_comm.recv().m_row.m_dst_onv.is_added_to_row());
    }

//    void on_row_send_(size_t irow) override {
//
//    }
//
//    void on_row_recv_(size_t irow) override {
//
//    }

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
            const defs::wf_t& weight = row.m_weight(m_ipart);
            res += std::pow(weight, 2.0);
        }
        return mpi::all_sum(res);
    }

    defs::wf_comp_t l1_norm() const {
        defs::wf_comp_t res = 0.0;
        auto& row = m_store.m_row;
        for (row.restart(); row.in_range(); row.step()) {
            const defs::wf_t& weight = row.m_weight(m_ipart);
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
        defs::wf_t& weight = row.m_weight(m_ipart);
        m_delta_nwalker(0, 0) += std::abs(new_weight);
        m_delta_nwalker(0, 0) -= std::abs(weight);
        m_delta_l2_norm_square(0, 0) += std::pow(std::abs(new_weight), 2.0);
        m_delta_l2_norm_square(0, 0) -= std::pow(std::abs(weight), 2.0);
        weight = new_weight;

        if (std::abs(new_weight) >= m_opts.nadd_initiator) grant_initiator_status();
        else revoke_initiator_status();
    }

    void change_weight(const defs::wf_t &delta) {
        set_weight(m_store.m_row.m_weight(m_ipart) + delta);
    }

    void scale_weight(const double &factor) {
        set_weight(factor*m_store.m_row.m_weight(m_ipart));
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

    // TODO: return a pair?
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
        row.m_src_initiator.put(initiator);
        row.m_src_deterministic.put(deterministic);
        row.m_dst_ipart = dst_ipart;
        return irow;
    }
};

#endif //M7_WAVEFUNCTION_H
