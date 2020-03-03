//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H


#include "src/hamiltonian/Hamiltonian.h"
#include "WalkerCommunicator.h"
#include "RankAllocator.h"
#include <iomanip>
#include <iostream>
#include <src/data/NumericView.h>
#include <src/io/FciqmcStatsFile.h>

class Propagator {
public:
    const InputOptions &m_input;
    const std::unique_ptr<Hamiltonian> &m_ham;
    const RankAllocator<Determinant> &m_rank_allocator;
    double m_tau;
    defs::ham_comp_t m_shift;
    bool vary_shift = false;
    defs::ham_comp_t m_largest_spawn_magnitude = 0;

    Propagator(const InputOptions &input,
               const std::unique_ptr<Hamiltonian> &ham,
               const RankAllocator<Determinant> &rank_allocator);

    void diagonal(const NumericView<defs::ham_comp_t> &hdiag, NumericView<defs::ham_t> &weight,
                  defs::ham_comp_t &delta_square_norm) const;

    void add_to_spawn_list(const Determinant &determinant, const defs::ham_t &weight,
                           bool flag_parent_initiator, TableArray<SpawnList> &spawn_list) const;

    virtual void off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                              const NumericView<bool> flag_deterministic,
                              const NumericView<bool> flag_initiator, TableArray<SpawnList> &spawn_list) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    void update(const size_t icycle, defs::ham_comp_t norm, defs::ham_comp_t norm_growth) {
        if (icycle % m_input.shift_update_period) return;
        if (!vary_shift) {
            if (norm < m_input.nwalker_target) return;
            else vary_shift = true;
        }
        m_shift -= m_input.shift_damp * consts::real_log(norm_growth) / m_tau;
    }

    void write_iter_stats(FciqmcStatsFile &stats_file) {
        stats_file.m_timestep->write(m_tau);
        stats_file.m_diagonal_shift->write(m_shift);
    }


//void evolve(const Perforable)
};


#endif //M7_PROPAGATOR_H
