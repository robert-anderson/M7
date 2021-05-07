//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include <src/core/io/InteractiveVariable.h>
#include <queue>
#include "src/core/hamiltonian/Hamiltonian.h"
#include "Wavefunction.h"
#include "MagnitudeLogger.h"

class Reweighter {
    const buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_const_shifts;
    std::vector<std::queue<defs::ham_t>> m_xc_queues;

    Reweighter(const fields::Numbers<defs::ham_comp_t, defs::ndim_wf>& const_shifts):
            m_const_shifts(const_shifts){}

};


struct Shift {
    const Options &m_opts;
    buffered::Numbers<defs::wf_t, defs::ndim_wf> m_nwalker_last_update;
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_values;
    Epochs m_variable_mode;
    InteractiveVariable<defs::wf_comp_t> m_nwalker_target;

    Shift(const Options &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    const defs::ham_comp_t & operator[](const size_t& ipart);

    void update(const Wavefunction& wf, const size_t& icycle, const double& tau);
};


class Propagator {
public:
    const NdFormat<defs::ndim_wf> &m_wf_fmt;
    const Hamiltonian<> &m_ham;
    const Options &m_opts;
    MagnitudeLogger m_magnitude_logger;
    Shift m_shift;
    /*
     * working objects
     */
    mutable buffered::Onv<> m_dst_onv;
    mutable conn::Antisym<> m_aconn;
    mutable OccupiedOrbitals m_occ;
    mutable VacantOrbitals m_vac;

    Propagator(const Options &opts, const Hamiltonian<> &ham, const NdFormat<defs::ndim_wf> &wf_fmt) :
            m_wf_fmt(wf_fmt),
            m_ham(ham),
            m_opts(opts),
            m_magnitude_logger(opts, ham.nsite(), ham.nelec()),
            m_shift(opts, wf_fmt),
            m_dst_onv(ham.nsite()),
            m_aconn(ham.nsite()),
            m_occ(ham.nsite()),
            m_vac(ham.nsite()) {}

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

};

#endif //M7_PROPAGATOR_H
