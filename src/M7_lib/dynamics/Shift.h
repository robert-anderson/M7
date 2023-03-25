//
// Created by Robert J. Anderson on 03/06/2021.
//

#ifndef M7_SHIFT_H
#define M7_SHIFT_H

#include <queue>

#include <M7_lib/io/InteractiveVariable.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/parallel/Epoch.h>

#include "MagnitudeLogger.h"
#include "M7_lib/wavefunction/Reference.h"
#include "M7_lib/wavefunction/Wavefunction.h"

/**
 * responsible for defining and updating the shift subtracted from the diagonal H elements in
 * propagation
 */
struct Shift {
    const conf::Document &m_opts;
    /**
     * the numbers of walkers on each WF part in the last period is stored so that the growth rate can be computed
     */
    buffered::Numbers<wf_comp_t, c_ndim_wf> m_nwalker_last_period;
    /**
     * values of the diagonal shift for each WF part
     */
    buffered::Numbers<ham_comp_t, c_ndim_wf> m_values;
    /**
     * queues storing the most recent values of the shift in order to enable constant time
     * computation of the rolling average shift values
     */
    v_t<std::queue<ham_comp_t>> m_avg_value_histories;
    /**
     * total unnormalized average shift values
     */
    buffered::Numbers<ham_comp_t, c_ndim_wf> m_avg_values;
    /**
     * the shift is initially held constant until the L1 norm of the wavefunction reaches
     * a user-defined level, at this point the shift is allowed to vary such that the L1
     * norm is stabilized
     */
    Epochs m_variable_mode;
    /**
     * the target walker number can be changed by the user on the fly
     */
    InteractiveVariable<wf_comp_t> m_nwalker_target;
    /**
     * if using target-driven damping, this will be y^2/4 where y is the normal (static) damp factor, else it will be 0
     */
    const double m_target_damp_fac;

    Shift(const conf::Document &opts, const NdFormat<c_ndim_wf>& wf_fmt);

    const ham_comp_t & operator[](uint_t ipart);

    /**
     * compute the change in all parts of the shift value based on the current values of wf.m_nwalkers
     *
     * the shift value is changed for a WF part if:
     *  - the walker number for the part exceeds the target
     *  - the cycle number is a nonzero multiple of the configured update period
     *
     * @param wf
     *  wavefunction whose population growth defines the change in shift
     * @param icycle
     *  MC cycle index
     * @param tau
     *  current timestep
     */
    void update(const wf::Vectors& wf, uint_t icycle, double tau, ham_comp_t value);

private:

    /**
     * add the current m_values array into the average queues and update the unnormalized averages accordingly
     */
    void add_to_average();

    /**
     * @param ipart
     *  part index
     * @return
     *  normalized average for ipart
     */
    ham_comp_t get_average(const uint_t& ipart) const;

};

#endif //M7_SHIFT_H
