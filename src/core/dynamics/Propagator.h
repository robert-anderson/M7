//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include <src/core/io/InteractiveVariable.h>
#include "src/core/hamiltonian/Hamiltonian.h"
#include "Wavefunction.h"
#include "MagnitudeLogger.h"


struct Shift {
    const Options &m_opts;
    buffered::Numbers<defs::wf_t, defs::ndim_wf> m_nwalker_last_update;
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_values;
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_avg_values; // un-normalised average (i.e. sum)
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_const_shift;
    Epochs m_variable_mode;
    Epochs m_reweighting_active;
    std::vector<std::queue<defs::ham_comp_t>> m_reweighting_factors;
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_total_reweighting;
    InteractiveVariable<defs::wf_comp_t> m_nwalker_target;

    Shift(const Options &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    const defs::ham_comp_t & operator[](const size_t& ipart);

    void update(const Wavefunction& wf, const size_t& icycle, const double& tau);

    void evaluate_reweighting(const size_t& npart, const size_t &icycle,
                                  const double& tau);
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
