//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>
#include "src/io/InputOptions.h"
#include "src/hamiltonian/Hamiltonian.h"
#include "Propagator.h"
#include "Wavefunction.h"

class FciqmcCalculation {
    const InputOptions m_input;
    Propagator m_prop;
    Wavefunction m_psi;
public:
/*
    FciqmcCalculation(const InputOptions &input) :
            m_input(input),
            m_h(input.fcidump_path),
            m_p(m_h),
            m_wf(input, m_h.choose_reference(input.spin_level)) {
        m_wf.m_reference.print();
    }*/

    void execute(size_t ncycle);

};


#endif //M7_FCIQMCCALCULATION_H
