//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_PROPAGATOR_H
#define M7_PROPAGATOR_H

#include "src/core/hamiltonian/Hamiltonian.h"
#include "Wavefunction.h"
#include "MagnitudeLogger.h"

class Propagator {
public:
    const Hamiltonian &m_ham;
    const Options &m_opts;
    MagnitudeLogger m_magnitude_logger;
    defs::ham_comp_t m_shift;
    /*
     * working objects
     */
    mutable elements::Onv m_dst_onv;
    mutable conn::AsOnv m_aconn;
    mutable OccupiedOrbitals m_occ;
    mutable VacantOrbitals m_vac;

//    Epoch& m_variable_shift;
//    Epoch& m_semi_stochastic;

    Propagator(const Hamiltonian &ham, const Options& opts):
    m_ham(ham),
    m_opts(opts),
    m_magnitude_logger(opts, ham.nsite(), ham.nelec()),
    m_shift(m_opts.shift_initial),
    m_dst_onv(ham.nsite()),
    m_aconn(m_dst_onv),
    m_occ(m_dst_onv),
    m_vac(m_dst_onv)
    {}


    //Propagator(FciqmcCalculation *fciqmc);

    //virtual ~Propagator()= default;

    virtual void diagonal(Wavefunction& m_wf, const size_t& irow) = 0;

    virtual void off_diagonal(Wavefunction& m_wf, const size_t& irow) = 0;

    virtual defs::ham_t round(const defs::ham_t &weight) {
        return weight;
    }

//    Epoch &variable_shift() {
//        return m_variable_shift;
//    }
//
//    Epoch &semi_stochastic() {
//        return m_semi_stochastic;
//    }

    const double& tau() const {
        return m_magnitude_logger.m_tau;
    }

//    void update(const size_t& icycle, defs::wf_comp_t nwalker, defs::wf_comp_t nwalker_growth);
//
//    void write_iter_stats(FciqmcStatsFile* stats_file) {
//        if (!mpi::i_am_root()) return;
//        stats_file->m_timestep.write(tau());
//        stats_file->m_diagonal_shift.write(m_shift);
//        stats_file->m_psingle.write(m_magnitude_logger.m_psingle);
//    }
};

#endif //M7_PROPAGATOR_H