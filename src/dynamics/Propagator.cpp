//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"


Propagator::Propagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                       const RankAllocator<Determinant> &rank_allocator) :
        m_input(input), m_ham(ham), m_rank_allocator(rank_allocator) {
    m_tau = input.tau_initial;
    m_shift = input.shift_initial;
}


void Propagator::diagonal(const NumericView<defs::ham_comp_t> &hdiag, NumericView<defs::ham_t> &weight,
                          defs::ham_comp_t &delta_square_norm) const {
    auto tmp = *weight;
    weight *= (1.0 - (hdiag - m_shift) * m_tau);
    delta_square_norm += std::pow(std::abs(*weight), 2) - std::pow(std::abs(tmp), 2);
}

void Propagator::add_to_spawn_list(const Determinant &determinant, const defs::ham_t &weight,
                                   bool flag_parent_initiator, TableArray<SpawnList> &spawn_lists) const {
    if (consts::float_is_zero(weight)) return;
    auto &spawn_list = spawn_lists[m_rank_allocator.get_rank(determinant)];
    spawn_list.push(determinant, weight, flag_parent_initiator);
}
