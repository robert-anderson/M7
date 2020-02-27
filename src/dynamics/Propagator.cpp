//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"


Propagator::Propagator(const Hamiltonian &h) : m_h(h) {}

void Propagator::diagonal(const NumericView<defs::ham_comp_t> &hdiag, NumericView<defs::ham_t> &weight,
                          defs::ham_comp_t &delta_square_norm) const {
    auto tmp = *weight;
    weight *= (1.0 - (hdiag - m_shift) * m_tau);
    delta_square_norm += std::pow(std::abs(weight), 2) - std::pow(std::abs(tmp), 2);
}
