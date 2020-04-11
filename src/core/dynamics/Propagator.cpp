//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"


Propagator::Propagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                       const RankAllocator<DeterminantElement> &rank_allocator) :
    m_input(input), m_ham(ham), m_rank_allocator(rank_allocator) {
    m_tau = input.tau_initial;
    m_shift = input.shift_initial;
}


void Propagator::diagonal(
    const NumericElement<defs::ham_comp_t> &hdiag,
    NumericElement<defs::ham_t> &weight, defs::ham_comp_t &delta_square_norm) const {
    auto tmp = weight;
    weight *= (1.0 - (*hdiag - m_shift) * m_tau);
    delta_square_norm += std::pow(std::abs(*weight), 2) - std::pow(std::abs(*tmp), 2);
}
