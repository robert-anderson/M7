//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include "src/core/hamiltonian/Hamiltonian.h"
#include "Wavefunction.h"
#include "MagnitudeLogger.h"

#if 0
class Propagator {
public:
    const Hamiltonian<> &m_ham;
    const Options &m_opts;
    MagnitudeLogger m_magnitude_logger;
    defs::ham_comp_t m_shift;
    /*
     * working objects
     */
    mutable elements::Onv<> m_dst_onv;
    mutable conn::Antisym<> m_aconn;
    mutable OccupiedOrbitals m_occ;
    mutable VacantOrbitals m_vac;

    Epoch m_variable_shift;

    Propagator(const Options &opts, const Hamiltonian<> &ham) :
    m_ham(ham),
    m_opts(opts),
    m_magnitude_logger(opts, ham.nsite(), ham.nelec()),
    m_shift(m_opts.shift_initial),
    m_dst_onv(ham.nsite()),
    m_aconn(m_dst_onv),
    m_occ(m_dst_onv),
    m_vac(m_dst_onv),
    m_variable_shift("variable shift mode"){}

    virtual void diagonal(Wavefunction& m_wf, const size_t& irow) = 0;

    virtual void off_diagonal(Wavefunction& m_wf, const size_t& irow) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

    const double& tau() const {
        return m_magnitude_logger.m_tau;
    }

    void update(const size_t& icycle, const Wavefunction& wf);

};

#endif //M7_PROPAGATOR_H
#endif //M7_PROPAGATOR_H