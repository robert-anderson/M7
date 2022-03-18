//
// Created by rja on 03/06/2021.
//

#ifndef M7_SHIFT_H
#define M7_SHIFT_H

#include <io/InteractiveVariable.h>
#include <queue>
#include <hamiltonian/Hamiltonian.h>
#include <wavefunction/Wavefunction.h>

#include "MagnitudeLogger.h"

/**
 * handles data and methods relating to the reweighting adaptation which corrects for
 * population control bias as described in Phys. Rev. B 103, 155135
 */
struct Reweighter {
    const fciqmc_config::Shift &m_opts;
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_const_shift;
    /**
     * the reweighting epochs begin some specified number of cycles after variable shift mode begins
     */
    Epochs m_active;
    /**
     * values of exp(tau* (const_shift - inst_shift)) for at most the past n cycles
     * where n is m_opts.ncycle_reweight_lookback.
     */
    std::vector<std::queue<defs::ham_comp_t>> m_histories;
    /**
     * products over all histories
     */
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_total;
    Reweighter(const fciqmc_config::Shift& opts, const NdFormat<defs::ndim_wf>& wf_fmt);
    /**
     * decide whether the reweighting epoch should begin, and if it should, fix the constant shift to the current
     * average shift value
     * @param icycle
     *  MC cycle index
     * @param ipart
     *  part index
     * @param begin_cond
     *  true if the epoch should begin if it is not already active
     * @param av_shift
     *  average shift value over some number of historical instantaneous values
     */
    void update(size_t icycle, size_t ipart, bool begin_cond, defs::ham_comp_t av_shift);
    /**
     * if the reweighting adapation is in use, and the epoch is active, then add to histories
     * @param ipart
     *  part index
     * @param shift
     *  instantaneous shift
     * @param icycle
     *  MC cycle index
     * @param tau
     *  current timestep
     */
    void add(size_t ipart, defs::ham_comp_t shift, size_t icycle, double tau);

private:
    /**
     * push new reweighting factor into the queue multiply it into the total. pop old factors and divide if the queue
     * size limit has been reached
     * @param ipart
     * @param v
     *  instantaneous shift
     */
    void add_to_history(size_t ipart, const defs::ham_comp_t& v);

    /**
     * @return
     *  true if the product over queues matches m_total (allowing for floating point non-commutivity)
     */
    bool product_chk() const;
};

/**
 * responsible for defining and updating the shift subtracted from the diagonal H elements in
 * propagation
 */
struct Shift {
    const fciqmc_config::Document &m_opts;
    /**
     * the numbers of walkers on each WF part in the last period is stored so that the growth rate can be computed
     */
    buffered::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_last_period;
    /**
     * values of the diagonal shift for each WF part
     */
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_values;
    /**
     * queues storing the most recent values of the shift in order to enable constant time
     * computation of the rolling average shift values
     */
    std::vector<std::queue<defs::ham_comp_t>> m_avg_value_histories;
    /**
     * total unnormalized average shift values
     */
    buffered::Numbers<defs::ham_comp_t, defs::ndim_wf> m_avg_values;
    /**
     * the shift is initially held constant until the L1 norm of the wavefunction reaches
     * a user-defined level, at this point the shift is allowed to vary such that the L1
     * norm is stabilized
     */
    Epochs m_variable_mode;
    /**
     * the target walker number can be changed by the user on the fly
     */
    InteractiveVariable<defs::wf_comp_t> m_nwalker_target;
    Reweighter m_reweighter;

    Shift(const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    const defs::ham_comp_t & operator[](const size_t& ipart);

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
    void update(const Wavefunction& wf, const size_t& icycle, const double& tau);

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
    defs::ham_comp_t get_average(const size_t& ipart) const;

};

#endif //M7_SHIFT_H
