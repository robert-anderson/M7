//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include <src/enumerators/VectorCombinationEnumerator.h>
#include "WalkerList.h"
#include "src/sample/HeatBathSampler.h"
#include "WalkerCommunicator.h"
#include "Propagator.h"
#include "RankAllocator.h"

class Wavefunction {
    WalkerList m_walker_list;
    WalkerCommunicator m_walker_communicator;

    Determinant m_reference;
    defs::ham_t m_reference_energy_numerator;
    defs::ham_t m_reference_energy_denominator;

public:

    defs::ham_comp_t m_square_norm;
    defs::ham_comp_t m_delta_square_norm;
    double m_norm_growth_rate = 0;

    Wavefunction(const std::unique_ptr<Propagator> &propagator, const Determinant &reference,
            double nwalker_initial, size_t nrow_walkers, size_t nrow_send, size_t nrow_recv);

    void propagate(const std::unique_ptr<Propagator> &propagator);

    void communicate();

    void consolidate_incoming_weight(){
        // TODO
    }

    void annihilate(const std::unique_ptr<Propagator> &propagator);

    void write_iter_stats(){
        std::cout << "=====   " << m_reference_energy_numerator/m_reference_energy_denominator << std::endl;
    }

    defs::ham_comp_t norm() const;

};


#endif //M7_WAVEFUNCTION_H
