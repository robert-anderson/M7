//
// Created by rja on 18/05/2020.
//

#include "MagnitudeLogger.h"

MagnitudeLogger::MagnitudeLogger(const Options &input) :
        m_input(input), m_tau(m_input.tau_initial) {}

void MagnitudeLogger::log(size_t nexcit, defs::ham_t helem, defs::prob_t prob) {
    defs::ham_comp_t tmp_hi_mag;
    if (nexcit == 1) {
        m_priv_nsingle.get()++;
        tmp_hi_mag = std::abs(helem) / prob;
        auto hi_mag = m_priv_hi_mag_single.get();
        if (tmp_hi_mag > hi_mag) hi_mag = tmp_hi_mag;
    } else if (nexcit == 2) {
        m_priv_ndouble.get()++;
        tmp_hi_mag = std::abs(helem) / prob;
        auto hi_mag = m_priv_hi_mag_double.get();
        if (tmp_hi_mag > hi_mag) hi_mag = tmp_hi_mag;
    }
}

void MagnitudeLogger::synchronize() {
    // thread reduce
    /*
    reduction::max(m_priv_hi_mag_single.get(), m_hi_mag_single);
    reduction::max(m_priv_hi_mag_double.get(), m_hi_mag_double);
    reduction::max(m_priv_nsingle.get(), m_nsingle);
    reduction::max(m_priv_ndouble.get(), m_ndouble);
     */
    m_hi_mag_single = m_priv_hi_mag_single.get();
    m_hi_mag_double = m_priv_hi_mag_double.get();
    m_nsingle = m_priv_nsingle.get();
    m_ndouble = m_priv_ndouble.get();

    // process reduce
    m_hi_mag_single = mpi::all_max(m_hi_mag_single);
    m_hi_mag_double = mpi::all_max(m_hi_mag_double);
    if (m_input.dynamic_tau){
        if (!m_enough_singles_for_dynamic_tau){
            m_enough_singles_for_dynamic_tau = m_nsingle > m_input.nenough_spawns_for_dynamic_tau;
            m_enough_singles_for_dynamic_tau = mpi::all_land(m_enough_singles_for_dynamic_tau);
        }
        if (!m_enough_doubles_for_dynamic_tau){
            m_enough_doubles_for_dynamic_tau = m_ndouble > m_input.nenough_spawns_for_dynamic_tau;
            m_enough_doubles_for_dynamic_tau = mpi::all_land(m_enough_doubles_for_dynamic_tau);
        }

        defs::ham_comp_t hi_mag_sum = m_hi_mag_single + m_hi_mag_double;
        ASSERT(hi_mag_sum>0.0);
        /*
         * highest transferred weight ~ tau x max(helem/prob)
         * i.e. recommended tau = max bloom / max(helem/prob)
         */
        m_tau = m_input.max_bloom/hi_mag_sum;
        m_psingle = m_hi_mag_single/hi_mag_sum;
    }
}
