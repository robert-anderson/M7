//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"


void Propagator::update(const size_t& icycle, const Wavefunction& wf) {
    //m_magnitude_logger.synchronize(icycle);
    m_shift.update(wf, icycle, tau());
}

void Propagator::load_fn(hdf5::GroupReader &parent) {
    size_t tmp;
    parent.load("nsite", tmp);
    parent.save("nelec", m_ham.nelec());
}

void Propagator::save_fn(hdf5::GroupWriter &parent) {
    parent.save("nsite", m_ham.nsite());
    parent.save("nelec", m_ham.nelec());
    parent.save("shift", m_magnitude_logger.m_tau);
    parent.save("tau", m_magnitude_logger.m_tau);
    parent.save("psingle", m_magnitude_logger.m_psingle);
}
