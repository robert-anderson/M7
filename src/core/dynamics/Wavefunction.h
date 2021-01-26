//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/io/Options.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include "src/core/table/CommunicatingPair.h"
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/field/Views.h"


struct Wavefunction : ra::Onv::Dynamic {
    const Options &m_opts;
    typedef BufferedTable<WalkerMappedTable> walkers_t;
    walkers_t m_walkers;
    typedef CommunicatingPair<SpawnTable> spawn_t;
    spawn_t m_spawn;

    ReductionSyndicate m_summables;

    ReductionMember<size_t, defs::ndim_wf> m_ninitiator;
    ReductionMember<int, defs::ndim_wf> m_delta_ninitiator;
    ReductionMember<size_t, defs::ndim_wf> m_nocc_onv;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_l2_norm_square;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_delta_l2_norm_square;

    Wavefunction(const Options &opts, size_t nsite, ra::Onv& ra) :
            ra::Onv::Dynamic(ra),
            m_opts(opts),
            // size_t nsite, size_t nroot, size_t nreplica, size_t nbucket
            m_walkers("walker table" , nsite, 1, 1, opts.nwalker_target / mpi::nrank()),
            // size_t nsite, size_t nroot, size_t nreplica
            m_spawn("spawning communicator", opts.spawn_buffer_expansion_factor, nsite, 1, 1),
            m_ninitiator(m_summables, {1, 1}),
            m_delta_ninitiator(m_summables, {1, 1}),
            m_nocc_onv(m_summables, {1, 1}),
            m_nwalker(m_summables, {1, 1}),
            m_delta_nwalker(m_summables, {1, 1}),
            m_l2_norm_square(m_summables, {1, 1}),
            m_delta_l2_norm_square(m_summables, {1, 1}) {
        m_walkers.set_expansion_factor(m_opts.spawn_buffer_expansion_factor);
    }

    void on_row_send_(size_t irow) override {

    }

    void on_row_recv_(size_t irow) override {

    }

    void begin_cycle() {
        m_summables.zero();
    }

    void end_cycle() {
        m_summables.all_sum();
    }

    void update(size_t icycle, double wait_time) {
        m_ra.update(icycle, wait_time, m_walkers, m_walkers.m_key_field);
    }

    void expand(size_t nrow_walker, size_t nrow_spawn) {
        m_walkers.expand(nrow_walker);
        m_spawn.expand(nrow_spawn);
    }

    defs::wf_comp_t square_norm() const {
        defs::wf_comp_t res = 0.0;
        for (size_t irow = 0; irow < m_walkers.m_hwm; ++irow) {
            if (!m_walkers.m_onv(irow).is_zero()) {
                res += std::pow(std::abs(m_walkers.m_weight(irow, 0, 0)), 2.0);
            }
        }
        return mpi::all_sum(res);
    }

    void grant_initiator_status(const size_t &irow) {
        auto view = m_walkers.m_flags.m_initiator(irow, 0, 0);
        if (!view) {
            m_delta_ninitiator(0, 0)++;
            view = true;
        }
    }

    void revoke_initiator_status(const size_t &irow) {
        auto view = m_walkers.m_flags.m_initiator(irow, 0, 0);
        if (view) {
            m_delta_ninitiator(0, 0)--;
            view = false;
        }
    }

    void set_weight(const size_t &irow, const defs::wf_t &new_weight) {
        m_delta_nwalker(0, 0) += std::abs(new_weight);
        m_delta_nwalker(0, 0) -= std::abs(m_walkers.m_weight(irow, 0, 0));
        m_delta_l2_norm_square(0, 0) += std::pow(std::abs(new_weight), 2.0);
        m_delta_l2_norm_square(0, 0) -= std::pow(std::abs(m_walkers.m_weight(irow, 0, 0)), 2.0);

        m_walkers.m_weight(irow, 0, 0) = new_weight;

        if (std::abs(new_weight) >= m_opts.nadd_initiator) grant_initiator_status(irow);
        else revoke_initiator_status(irow);
    }

    void change_weight(const size_t &irow, const defs::wf_t &delta) {
        set_weight(irow, m_walkers.m_weight(irow, 0, 0) + delta);
    }

    void scale_weight(const size_t &irow, const double &factor) {
        set_weight(irow, m_walkers.m_weight(irow, 0, 0) * factor);
    }

    void zero_weight(const size_t &irow) {
        set_weight(irow, 0.0);
    }

    void remove_walker(const size_t &irow) {
        const auto onv = m_walkers.m_onv(irow);
        auto lookup = m_walkers[onv];
        zero_weight(irow);
        m_walkers.erase(lookup);
    }

    size_t create_walker(const views::Onv<> &onv, const defs::ham_t weight,
                         const defs::ham_comp_t &hdiag, bool refconn) {
        ASSERT(mpi::i_am(m_ra.get_rank(onv)));
        if (m_walkers.is_full()) m_walkers.expand(1);
        auto irow = m_walkers.insert(onv);
        ASSERT(m_walkers.m_onv(irow) == onv)
        set_weight(irow, weight);
        m_walkers.m_hdiag(irow) = hdiag;
        m_walkers.m_flags.m_reference_connection(irow) = refconn;
        m_walkers.m_flags.m_deterministic(irow) = false;
        return irow;
    }

    // TODO: return a pair?
    size_t add_spawn(const views::Onv<> &dst_onv, const defs::wf_t &delta, bool initiator, bool deterministic) {
        auto irank = m_ra.get_rank(dst_onv);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "SENDING SPAWNED WALKER" << std::endl;
        std::cout << consts::verb << "generated determinant:   " << dst_onv.to_string() << std::endl;
        std::cout << consts::verb << "destination rank:        " << irank << std::endl;
        std::cout << consts::verb << "spawned weight:          " << delta << std::endl;
        std::cout << consts::verb << "parent is initiator:     " << initiator << std::endl;
        std::cout << consts::verb << "parent is deterministic: " << deterministic << std::endl;
#endif
        auto &dst_table = m_spawn.send(irank);

        auto irow = dst_table.push_back();
        dst_table.m_dst_onv(irow) = dst_onv;
        dst_table.m_delta_weight(irow, 0, 0) = delta;
        dst_table.m_flags.m_src_initiator(irow, 0, 0) = initiator;
        dst_table.m_flags.m_src_deterministic(irow) = deterministic;
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
