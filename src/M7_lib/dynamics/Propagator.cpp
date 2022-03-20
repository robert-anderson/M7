//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"


void Propagator::update(const size_t& icycle, const Wavefunction& wf) {
    m_shift.update(wf, icycle, tau());
}

void Propagator::load_fn(hdf5::GroupReader &parent) {
    REQUIRE_EQ_ALL(parent.load<size_t>("nsite"), m_ham.m_bd.m_nsite, "number of sites is not consistent with archive");
    REQUIRE_EQ_ALL(parent.load<size_t>("nmode"), m_ham.m_bd.m_nmode, "number of boson modes is not consistent with archive");
    REQUIRE_EQ_ALL(parent.load<size_t>("nelec"), m_ham.nelec(), "number of electrons is not consistent with archive");
}

void Propagator::save_fn(hdf5::GroupWriter &parent) {
    hdf5::GroupWriter gw("propagator", parent);
    gw.save("nsite", m_ham.m_bd.m_nsite);
    gw.save("nmode", m_ham.m_bd.m_nmode);
    gw.save("nelec", m_ham.nelec());
    //gw.save("shift", m_shift.m_values);
    gw.save("tau", m_tau);
    //gw.save("psingle", m_magnitude_logger.m_psingle);
}

void Propagator::imp_samp_delta(wf_t &delta, const Mbf &src_mbf, const Mbf &dst_mbf, const ham_comp_t &src_energy) const {
    auto& g = m_opts.m_propagator.m_imp_samp_exp.get();
    // check that importance sampling is in use at all before evaluating energies
    if (g) delta*=std::exp(g*(src_energy-m_ham.get_energy(dst_mbf)));
}

void Propagator::imp_samp_delta(wf_t &delta, const Mbf &src_mbf, const Mbf &dst_mbf) const {
    imp_samp_delta(delta, src_mbf, dst_mbf, m_ham.get_energy(src_mbf));
}
