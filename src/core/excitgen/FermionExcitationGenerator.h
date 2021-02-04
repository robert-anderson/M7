//
// Created by RJA on 20/11/2020.
//

#ifndef M7_FERMIONEXCITATIONGENERATOR_H
#define M7_FERMIONEXCITATIONGENERATOR_H

#include "ExcitationGenerator.h"

class FermionExcitationGenerator : public ExcitationGenerator {
protected:
    const size_t m_nexcit;
    const bool m_spin_conserving;

public:
    FermionExcitationGenerator(const Hamiltonian<> *h, PRNG &prng, size_t nexcit);
    
};


#endif //M7_FERMIONEXCITATIONGENERATOR_H
