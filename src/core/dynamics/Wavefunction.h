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
#include "src/core/field/Views.h"


struct Wavefunction {
    const Options &m_opts;
    typedef BufferedTable<WalkerTable> table_t;
    table_t m_walkers;
    typedef CommunicatingPair<SpawnTable> spawn_t;
    spawn_t m_spawn;
    typedef RankAllocator<fields::Onv> rank_alloc_t;
    rank_alloc_t m_ra;

    Wavefunction(const Options &opts, fields::Onv::params_t onv_params) :
            m_opts(opts),
            m_walkers("walker table", 1000, onv_params, 1, 1),
            m_spawn("spawning communicator", onv_params, 1, 1),
            m_ra(100, 10) {
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

    defs::wf_comp_t unnorm_energy(const Hamiltonian &ham) const {
        defs::wf_comp_t res = 0.0;
        for (size_t irow = 0; irow < m_walkers.m_hwm; ++irow) {
            if (m_walkers.m_onv(irow).is_zero()) continue;
            defs::wf_t weighti = consts::conj(m_walkers.m_weight(irow, 0, 0));
            for (size_t jrow = 0; jrow < m_walkers.m_hwm; ++jrow) {
                if (m_walkers.m_onv(jrow).is_zero()) continue;
                defs::wf_t weightj = m_walkers.m_weight(jrow, 0, 0);
                res += weighti * ham.get_element(m_walkers.m_onv(irow), m_walkers.m_onv(jrow)) * weightj;
            }
        }
        return mpi::all_sum(res);
    }

    defs::wf_comp_t energy(const Hamiltonian &ham) const {
        return unnorm_energy(ham)/square_norm();
    }

    size_t add_walker(const views::Onv &onv, const defs::ham_t weight, const defs::ham_comp_t &hdiag,
                      bool refconn, bool initiator) {
        auto irow = m_walkers.insert(onv);
        ASSERT(m_walkers.m_onv(irow) == onv)
        m_walkers.m_weight(irow, 0, 0) = weight;
        m_walkers.m_hdiag(irow) = hdiag;
        m_walkers.m_flags.m_reference_connection(irow) = refconn;
        m_walkers.m_flags.m_deterministic(irow) = false;
        m_walkers.m_flags.m_initiator(irow, 0, 0) = initiator;
        return irow;
    }

    // TODO: return a pair?
    size_t add_spawn(const views::Onv &dst_onv, const defs::wf_t &delta, bool initiator, bool deterministic) {
        auto irank = m_ra.get_rank(dst_onv);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "SENDING SPAWNED WALKER" << std::endl;
        std::cout << consts::verb << "generated determinant:   " << dst_onv.to_string() << std::endl;
        std::cout << consts::verb << "destination rank:        " << irank << std::endl;
        std::cout << consts::verb << "spawned weight:          " << delta << std::endl;
        std::cout << consts::verb << "parent is initiator:     " << flag_initiator << std::endl;
        std::cout << consts::verb << "parent is deterministic: " << flag_deterministic << std::endl;
#endif
        auto &dst_table = m_spawn.send(irank);
        if (dst_table.m_hwm + 1 == dst_table.m_nrow) {
            /*
             * the table is full
             */
            std::cout << "Spawn send table is full, reallocating..." << std::endl;
            m_spawn.send().expand(0.5 * dst_table.m_nrow);
        }
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
