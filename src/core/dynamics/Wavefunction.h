//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/io/ParallelizationStatsFile.h>
#include "src/core/fermion/Determinant.h"
#include "src/core/io/Options.h"
#include "src/core/util/Timer.h"
#include "WalkerList.h"
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

    WalkerList m_data;
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


private:

    defs::ham_comp_t projected_energy(Hamiltonian *ham, const DeterminantElement &ref) {
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

    defs::wf_comp_t nwalker() {
        // debug only
        Reducible<defs::wf_comp_t> nw;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            if (m_data.row_empty(i)) continue;
            nw += std::abs(*m_data.m_weight(i));
        }
        nw.mpi_sum();
        return nw.reduced();
    }

    size_t nocc() {
        // debug only
        Reducible<size_t> nocc;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            nocc+=(size_t)!m_data.row_empty(i);
        }
        nocc.mpi_sum();
        return nocc.reduced();
    }

    size_t ninitiator() {
        // debug only
        Reducible<size_t> ninit;
        for (size_t i = 0ul; i < m_data.high_water_mark(0); ++i) {
            if (m_data.row_empty(i)) continue;
            ninit+=(size_t)m_data.m_flags.m_initiator(i);
        }
        ninit.mpi_sum();
        return ninit.reduced();
    }

    void annihilate_row(const size_t &irow_recv);

};


#endif //M7_WAVEFUNCTION_H
