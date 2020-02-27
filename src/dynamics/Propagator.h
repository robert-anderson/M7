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

class Propagator {
public:
    const std::unique_ptr<Hamiltonian> &m_ham;
    const RankAllocator<Determinant> &m_rank_allocator;
    double m_tau;
    defs::ham_comp_t m_shift;
    defs::ham_comp_t m_largest_spawn_magnitude;

    Propagator(const std::unique_ptr<Hamiltonian> &ham,
               const RankAllocator<Determinant> &rank_allocator, double tau, defs::ham_comp_t shift);

    void diagonal(const NumericView<defs::ham_comp_t> &hdiag, NumericView<defs::ham_t> &weight,
                  defs::ham_comp_t &delta_square_norm) const;

    void add_to_spawn_list(const Determinant &determinant, const defs::ham_t &weight,
                           bool flag_parent_initiator, TableArray<SpawnList> &spawn_list) const;

    virtual void off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                              const NumericView<bool> flag_deterministic,
                              const NumericView<bool> flag_initiator, TableArray<SpawnList> &spawn_list) const = 0;

    void update(const size_t icycle, defs::ham_comp_t norm_growth) {
        m_shift -= consts::real_log(norm_growth) / m_tau;
        std::cout << "shift: " << std::setprecision(10) << m_shift << std::endl;
    }

    //void evolve(const Perforable)
};


#endif //M7_PROPAGATOR_H
