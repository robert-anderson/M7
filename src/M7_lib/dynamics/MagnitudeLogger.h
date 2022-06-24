//
// Created by Robert J. Anderson on 18/05/2020.
//

#ifndef M7_MAGNITUDELOGGER_H
#define M7_MAGNITUDELOGGER_H


#include <cstddef>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/Reduction.h>
#include <M7_lib/excitgen/ExcitGenGroup.h>
#include <M7_lib/conf/Conf.h>

struct MagnitudeLogger {
    /**
     * desired maximum number of walkers spawned from an integer walker
     */
    const ham_comp_t m_max_bloom;
    /**
     * number of draws required to be made from each case before it is allowed to affect dynamic probabilities and tau
     */
    const uint_t m_ndraw_min;
    const bool m_static_tau;
    const bool m_static_probs;
    const double m_tau_min, m_tau_max, m_prob_min;
    const uint_t m_period;
    /**
     * number of draws made in each excitation case
     */
    NdReduction<ham_comp_t, 1ul> m_ndraw;
    /**
     * largest generated magnitude (|helement|/prob) from a single integerized walker
     */
    NdReduction<ham_comp_t, 1ul> m_gamma;
    /**
     * working array for assigning new probability values to excitation generator group object
     */
    mutable std::vector<prob_t> m_new_probs;
    MagnitudeLogger(ham_comp_t max_bloom, uint_t ndraw_min, uint_t nexcase, bool static_tau, bool static_probs,
                    double tau_min, double tau_max, double prob_min, uint_t period);

    /**
     * add a generated spawn to the logger, i.e. increment the draw counter, and update the highest magnitude if the
     * given event has a higher magnitude than the current value
     * @param icase
     *  index of the excitation level generated
     * @param helem
     *  off-diagonal Hamiltonian matrix element associated with the drawn excitation
     * @param prob
     *  probability that this excitation was drawn *given* icase (i.e. prob of attempting an icase excitation is not
     *  factored in)
     */
    void log(uint_t icase, const ham_t& helem, const prob_t& prob);

private:
    /**
     * excitation level x:
     *  max_mag_x = maximum spawned magnitude for x = tau * gamma_x / P_x                                          (*)
     * where gamma_x is the maximum value of helem/prob and P_x is the probability with which x was chosen
     *
     * updating the probabilities:
     * the aim is to even out max_mags i.e. max_mag_1 = max_mag_2 = ... = max_mag_n i.e. n-1 constraints
     * subject to the additional constraint that the sum of P_x is 1.0
     *
     * solving this linear system results in the following new value for the probabilities:
     *  P_x = gamma_x / (gamma_1 + gamma_2 + ... + gamma_n)
     *
     * choosing this value for each probability also makes the system of equations (*) degenerate, and so given an ideal
     * value for the maximum spawned magnitude (bloom), the timestep can simply be updated as:
     *  tau = max_bloom / (gamma_1 + gamma_2 + ... + gamma_n)
     *
     * @param tau
     *  reference to timestep value to update
     */
    void update_tau(double& tau, const ham_comp_t& gamma_sum);

public:

    void update(uint_t icycle, double& tau);

    /**
     * @param tau
     *  reference to timestep value to update
     * @param excit_gens
     *  reference to group of excitation generators to update probabilities for
     */
    void update(uint_t icycle, double& tau, ExcitGenGroup& excit_gens);
};

#endif //M7_MAGNITUDELOGGER_H
