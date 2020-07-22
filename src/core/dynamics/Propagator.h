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
#include <src/core/parallel/Reducible.h>
#include <src/core/parallel/Epoch.h>

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
    Reducible<defs::wf_comp_t> m_largest_spawn_magnitude;

    mutable Determinant m_dst_det;
    mutable AntisymConnection m_aconn;
    mutable OccupiedOrbitals m_occ;
    mutable VacantOrbitals m_vac;

    Epoch& m_variable_shift;
    Epoch& m_semi_stochastic;

    Propagator(FciqmcCalculation *fciqmc);

    virtual ~Propagator()= default;

    void spawn(SpawnList &spawn_list, const DeterminantElement &dst_det, const defs::wf_t &delta,
            bool flag_initiator, bool flag_deterministic) {
        auto irank = m_rank_allocator.get_rank(dst_det);
        spawn_list.add(irank, dst_det, delta, flag_initiator, flag_deterministic);
    }

    virtual void diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                          bool flag_deterministic,
                          defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) = 0;


    virtual void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                              SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    Epoch &variable_shift() {
        return m_variable_shift;
    }

    Epoch &semi_stochastic() {
        return m_semi_stochastic;
    }

    void update(const size_t icycle, defs::wf_comp_t nwalker, defs::wf_comp_t nwalker_growth) {
        m_magnitude_logger.synchronize();
        // TODO
        m_largest_spawn_magnitude.mpi_max();
        if (icycle % m_input.shift_update_period) return;
        m_variable_shift.update(icycle, nwalker >= m_input.nwalker_target);
        if (m_variable_shift) m_shift -= m_input.shift_damp * consts::real_log(nwalker_growth) / m_tau;
    }

    void write_iter_stats(FciqmcStatsFile* stats_file) {
        if (!mpi::i_am_root()) return;
        stats_file->m_timestep.write(m_tau);
        stats_file->m_diagonal_shift.write(m_shift);
        stats_file->m_psingle.write(m_magnitude_logger.m_psingle);
    }
};


#endif //M7_PROPAGATOR_H
