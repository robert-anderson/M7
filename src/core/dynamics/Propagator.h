//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H


#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/parallel/RankAllocator.h"
#include "SpawnList.h"
#include "MagnitudeLogger.h"
#include <iomanip>
#include <iostream>
#include <src/core/io/FciqmcStatsFile.h>
#include <src/core/parallel/Hybrid.h>

class FciqmcCalculation;

class Propagator {
public:
    FciqmcCalculation *m_fciqmc;
    const Options &m_input;
    const std::unique_ptr<Hamiltonian> &m_ham;
    const RankAllocator<DeterminantElement> &m_rank_allocator;
    MagnitudeLogger m_magnitude_logger;
    double m_tau;
    defs::ham_comp_t m_shift;
    Hybrid<defs::wf_comp_t> m_largest_spawn_magnitude;
    Distributed<size_t> m_icycle_vary_shift;

    Propagator(FciqmcCalculation *fciqmc);

    virtual ~Propagator()= default;

    void spawn(SpawnList &spawn_list, const DeterminantElement &dst_det, const defs::wf_t &delta, bool flag_initiator) {
        auto const mag = std::abs(delta);
        auto irank = m_rank_allocator.get_rank(dst_det);
        spawn_list.add(irank, dst_det, delta, flag_initiator);
    }

    virtual void diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                          bool flag_deterministic,
                          defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) = 0;


    virtual void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                              SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    bool varying_shift(){
        return m_icycle_vary_shift.reduced() != ~0ul;
    }

    const size_t &icycle_vary_shift(){
        return m_icycle_vary_shift.reduced();
    }

    void update(const size_t icycle, defs::wf_comp_t nwalker, defs::wf_comp_t nwalker_growth) {
        m_magnitude_logger.synchronize();
        m_largest_spawn_magnitude.max();
        if (icycle % m_input.shift_update_period) return;
        if (!varying_shift()) {
            if (nwalker >= m_input.nwalker_target) m_icycle_vary_shift = icycle;
            m_icycle_vary_shift.mpi_min();
        }
        if (varying_shift())
            m_shift -= m_input.shift_damp * consts::real_log(nwalker_growth) / m_tau;
    }

    void write_iter_stats(FciqmcStatsFile* stats_file) {
        if (!mpi::i_am_root()) return;
        stats_file->m_timestep.write(m_tau);
        stats_file->m_diagonal_shift.write(m_shift);
        stats_file->m_psingle.write(m_magnitude_logger.m_psingle);
    }
};


#endif //M7_PROPAGATOR_H
