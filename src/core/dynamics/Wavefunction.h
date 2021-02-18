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
#include "src/core/fieldz/FieldsZ.h"

struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow> {

    const Options &m_opts;

    NdFormat<defs::ndim_wf> m_format;
    // current "part" i.e. the flat element of the format
    size_t m_ipart = 0ul;

    ReductionSyndicate m_summables;

    ReductionMember<size_t, defs::ndim_wf> m_ninitiator;
    ReductionMember<int, defs::ndim_wf> m_delta_ninitiator;
    ReductionMember<size_t, defs::ndim_wf> m_nocc_onv;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_l2_norm_square;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_l2_norm_square;

    Wavefunction(const Options &opts, size_t nsite) :
            Communicator<WalkerTableRow, SpawnTableRow>(
                    "walker",
                    opts.walker_buffer_expansion_factor,
                    opts.nload_balance_block_per_rank*mpi::nrank(),
                    opts.load_balance_period,
                    WalkerTableRow(nsite, opts.nroot, opts.nreplica), SpawnTableRow(nsite),
                    MappedTableBaseZ::nbucket_guess(opts.nwalker_target/mpi::nrank(), 3),
                    opts.acceptable_load_imbalance
            ),
            m_opts(opts),
            m_format({opts.nroot, opts.nreplica}),
            m_ninitiator(m_summables, m_format),
            m_delta_ninitiator(m_summables, m_format),
            m_nocc_onv(m_summables, m_format),
            m_nwalker(m_summables, m_format),
            m_delta_nwalker(m_summables, m_format),
            m_l2_norm_square(m_summables, m_format),
            m_delta_l2_norm_square(m_summables, m_format) {
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
        const defs::wf_t& weight = row.m_weight(m_ipart);
        for (row.restart(); row.in_range(); row.step()) {
            res += std::pow(weight, 2.0);
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
        zero_weight();
        m_store.erase(lookup);
    }

    size_t create_walker_(const fieldsz::Onv<> &onv, const defs::ham_t weight,
                          const defs::ham_comp_t &hdiag, bool refconn) {
        ASSERT(mpi::i_am(m_ra.get_rank(onv)));
        if (m_store.is_full()) m_store.expand(1);
        auto irow = m_store.insert(onv);
        m_store.m_row.jump(irow);
        ASSERT(m_store.m_row.m_onv == onv)
        set_weight(weight);
        m_store.m_row.m_hdiag = hdiag;
        m_store.m_row.m_reference_connection.put(m_ipart, refconn);
        m_store.m_row.m_deterministic.clr(m_ipart);
        return irow;
    }


    Table::Loc create_walker(const fieldsz::Onv<> &onv, const defs::ham_t weight,
                          const defs::ham_comp_t &hdiag, bool refconn) {
        size_t irank = m_ra.get_rank(onv);
        size_t irow;
        if (mpi::i_am(irank)) irow = create_walker_(onv, weight, hdiag, refconn);
        mpi::bcast(irow, irank);
        return {irank, irow};
    }

    // TODO: return a pair?
    size_t add_spawn(const fieldsz::Onv<> &dst_onv, const defs::wf_t &delta,
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

#if 0

#include "src/core/io/ParallelizationStatsFile.h"
#include "src/core/observables/KramersSectorOccupation.h"
#include "src/core/basis/FermionOnv.h"
#include "src/core/io/Options.h"
#include "src/core/util/Timer.h"
#include "WalkerTable.h"
#include "SpawnList.h"
#include "Propagator.h"
#include "src/core/parallel/DistributedAccumulation.h"
#include "DeterministicSubspace.h"
#include "Reference.h"

class FciqmcCalculation;

class Wavefunction {
    FciqmcCalculation *m_fciqmc;
    const Options &m_input;
    const std::unique_ptr<Propagator> &m_prop;
    std::unique_ptr<DeterministicSubspace> m_detsub = nullptr;
    std::unique_ptr<KramersSectorOccupation> m_mk_sums = nullptr;

    WalkerTable m_data;
    SpawnList m_send, m_recv;
    Reference m_reference;

    Reducible<defs::wf_t> m_aborted_weight;
    DistributedAccumulation<size_t, int64_t> m_ninitiator;

    /*
     * timers for the three main sections of the wavefunction evolution
     */
    Timer m_propagation_timer;
    Timer m_communication_timer;
    Timer m_annihilation_timer;

    size_t nrow_free = 0ul;

    ParallelizationStatsFile *parallel_stats_file();;

public:

    /*
     * Square norm is sum_i(|w_i|^2)
     */
    DistributedAccumulation<defs::wf_comp_t> m_square_norm;
    /*
     * Walker number is sum_i(|w_i|)
     */
    DistributedAccumulation<defs::wf_comp_t> m_nwalker;
    defs::wf_comp_t m_nwalker_growth_rate;

    /*
     * number of occupied determinants, the principal variable is distributed, since
     * rank-resolved data could inform load balancing
     */
    DistributedAccumulation<size_t, int64_t> m_nocc_det;

    explicit Wavefunction(FciqmcCalculation *fciqmc);

    ~Wavefunction() {
        std::cout << "# initiators: " << m_data.verify_ninitiator(m_input.nadd_initiator) << std::endl;
        m_data.report_top_weighted();
    }

    /**
     * Effect arithmetic updates on member variables (no communication)
     */
    void update(const size_t &icycle);

    void propagate();

    void communicate();

    void consolidate_incoming_weight() {
        // TODO
    }

    void annihilate();

    /**
     * Perform the necessary thread and MPI reductions on members
     */
    void synchronize();

    void write_iter_stats(FciqmcStatsFile *stats_file);

    defs::ham_comp_t ref_proj_energy(){
        return m_reference.proj_energy();
    }


private:
    void annihilate_row(const size_t &irow_recv);

    defs::ham_comp_t projected_energy_check(FermionHamiltonian *ham, const DeterminantElement &ref) {
        // debug only
        Reducible<defs::ham_t> e;
        Reducible<defs::wf_t> norm;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            if (m_data.row_empty(i)) continue;
            auto det = m_data.m_determinant(i);
            if (det == ref) norm = *m_data.m_weight(i);
            e += *m_data.m_weight(i) * ham->get_element(ref, m_data.m_determinant(i));
        }
        e.mpi_sum();
        norm.mpi_sum();
        return consts::real(e.reduced() / norm.reduced());
    }

    defs::wf_comp_t nwalker_check() {
        // debug only
        Reducible<defs::wf_comp_t> nw;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            if (m_data.row_empty(i)) continue;
            nw += std::abs(*m_data.m_weight(i));
        }
        nw.mpi_sum();
        return nw.reduced();
    }

    size_t nocc_check() {
        // debug only
        Reducible<size_t> nocc;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            nocc+=(size_t)!m_data.row_empty(i);
        }
        nocc.mpi_sum();
        return nocc.reduced();
    }

    size_t ninitiator_check() {
        // debug only
        Reducible<size_t> ninit;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            if (m_data.row_empty(i)) continue;
            ninit+=(size_t)m_data.m_flags.m_initiator(i);
        }
        ninit.mpi_sum();
        return ninit.reduced();
    }

};


#endif //M7_WAVEFUNCTION_H
#endif //M7_WAVEFUNCTION_H
