//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"


void Propagator::update(const size_t& icycle, const Wavefunction& wf) {
    //m_magnitude_logger.synchronize(icycle);
    std::cout << bool(m_shift[0]) << std::endl;
    m_shift.update(wf, icycle, tau());
}