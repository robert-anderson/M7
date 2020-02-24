//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H


#include <src/integrals/AbInitioHamiltonian.h>
#include <iomanip>

class Propagator {
public:
    const AbInitioHamiltonian &m_h;
    double tau = 0.01;
    defs::ham_t m_shift = 0.0;
    Propagator(const AbInitioHamiltonian &h);

    void update_shift(defs::ham_comp_t norm_growth){
        m_shift -= consts::real_log(norm_growth)/tau;
        std::cout << "shift: " << std::setprecision(10) << m_shift <<std::endl;
    }

    //void evolve(const Perforable)
};


#endif //M7_PROPAGATOR_H
