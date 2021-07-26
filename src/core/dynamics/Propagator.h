//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include <src/core/config/FciqmcConfig.h>
#include <src/core/io/Archivable.h>
#include <src/core/basis/Suites.h>
#include "Shift.h"

/**
 * This class is not optionally archivable,
 * if an archive is being saved or checkpointed, the propagator data must be included
 * if an archive is being loaded, the propagator data must be read and verified against the data on file to ensure
 * compatibility with any objects reinstated from disk
 */
class Propagator : public Archivable {
public:
    const NdFormat<defs::ndim_wf> &m_wf_fmt;
    const Hamiltonian<> &m_ham;
    const fciqmc_config::Document &m_opts;
    MagnitudeLogger m_magnitude_logger;
    Shift m_shift;
    /*
     * working objects
     */
    mutable suite::Mbfs m_dst;
    mutable suite::Conns m_conn;
    mutable OccupiedOrbitals m_occ;
    mutable VacantOrbitals m_vac;

    Propagator(const fciqmc_config::Document &opts, const Hamiltonian<> &ham, const NdFormat<defs::ndim_wf> &wf_fmt) :
            Archivable("propagator", opts.m_archive),
            m_wf_fmt(wf_fmt),
            m_ham(ham),
            m_opts(opts),
            m_magnitude_logger(opts.m_propagator, ham.nsite(), ham.nelec()),
            m_shift(opts, wf_fmt),
            m_dst(ham.nsite()),
            m_conn(ham.nsite()),
            m_occ(ham.nsite()),
            m_vac(ham.nsite()) {}

    virtual ~Propagator() {}

    virtual void diagonal(Wavefunction &wf, const size_t &ipart) = 0;

    virtual void off_diagonal(Wavefunction &wf, const size_t &ipart) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    const double &tau() const {
        return m_magnitude_logger.m_tau;
    }

    void update(const size_t &icycle, const Wavefunction &wf);

    virtual bool is_exact() const = 0;

private:
    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;

};

#endif //M7_PROPAGATOR_H
