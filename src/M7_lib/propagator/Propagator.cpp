//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"

std::unique_ptr<guide::Wavefunction> Propagator::make_imp_samp_guide(const conf::GuideWavefunction& opts) const {
    const str_t fmt = "Guiding propagation via importance sampling with {} wavefunction";
    if (opts.m_gutzwiller_like.m_enabled) {
        logging::info(fmt, "Gutzwiller-like");
        return smart_ptr::make_poly_unique<guide::Wavefunction, guide::GutzwillerLike>(
                m_ham, opts.m_gutzwiller_like.m_fac.m_value);
    }
    else if (opts.m_suppress_multi_occ.m_enabled) {
        logging::info(fmt, "multiple occupation suppressing");
        return smart_ptr::make_poly_unique<guide::Wavefunction, guide::SuppressMultiOcc>(
                opts.m_suppress_multi_occ.m_fac.m_value);
    }
    return nullptr;
}

void Propagator::update(const uint_t& icycle, const Wavefunction& wf) {
    m_shift.update(wf, icycle, tau());
}

void Propagator::load_fn(const hdf5::NodeReader& parent) {
    REQUIRE_EQ_ALL(parent.load<uint_t>("nsite"), m_ham.m_frm.m_basis.m_nsite,
                   "number of fermion sites is not consistent with HDF5 archive");
    REQUIRE_EQ_ALL(parent.load<uint_t>("nmode"), m_ham.m_bos.m_basis.m_nmode,
                   "number of boson modes is not consistent with HDF5 archive");
}

void Propagator::save_fn(const hdf5::NodeWriter& parent) {
    hdf5::GroupWriter gw(parent, "propagator");
    gw.save("nsite", uint_t(m_ham.m_frm.m_basis.m_nsite));
    gw.save("nmode", m_ham.m_bos.m_basis.m_nmode);
    //gw.save("shift", m_shift.m_values);
    gw.save("tau", m_tau);
    //gw.save("psingle", m_magnitude_logger.m_psingle);
}

void Propagator::imp_samp_delta(wf_t& delta, const Mbf& dst_mbf, ham_comp_t src_energy) const {
    // check that importance sampling is in use at all before evaluating energies
    // if (m_imp_samp_guide) delta*=std::exp(m_imp_samp_exp*(src_energy-m_ham.get_energy(dst_mbf)));
}

void Propagator::imp_samp_delta(wf_t& delta, const Mbf& src_mbf, const Mbf& dst_mbf) const {
    if (m_imp_samp_guide) {
        delta*=m_imp_samp_guide->overlap(src_mbf);
        delta/=m_imp_samp_guide->overlap(dst_mbf);
    }
}
