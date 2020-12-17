//
// Created by rja on 18/05/2020.
//

#include "MagnitudeLogger.h"

MagnitudeLogger::MagnitudeLogger(const Options &input, defs::prob_t psingle) :
        m_input(input),
        m_enough_singles_for_dynamic_tau("use singles ratio for dynamic tau"),
        m_enough_doubles_for_dynamic_tau("use doubles ratio for dynamic tau"),
        m_psingle(psingle), m_tau(m_input.tau_initial) {}

void MagnitudeLogger::log(size_t nexcit, defs::ham_t helem, defs::prob_t prob) {
    defs::ham_comp_t tmp_hi_mag;
    if (nexcit == 1) {
        ++m_nsingle;
        tmp_hi_mag = std::abs(helem) / prob;
#ifdef VERBOSE_DEBUGGING
        if (m_tau*tmp_hi_mag>m_input.max_bloom) {
            std::cout << consts::verb << consts::chevs << "MAX BLOOM EXCEEDED BY SINGLE EXCITATION" << std::endl;
            std::cout << consts::verb << "weight transferred:     " << m_tau*tmp_hi_mag << std::endl;
        }
#endif
        if (tmp_hi_mag>m_hi_mag_single.local()) m_hi_mag_single = tmp_hi_mag;
    } else if (nexcit == 2) {
        ++m_ndouble;
        tmp_hi_mag = std::abs(helem) / prob;
#ifdef VERBOSE_DEBUGGING
        if (m_tau*tmp_hi_mag>m_input.max_bloom) {
            std::cout << consts::verb << consts::chevs << "MAX BLOOM EXCEEDED BY DOUBLE EXCITATION" << std::endl;
            std::cout << consts::verb << "weight transferred:     " << m_tau*tmp_hi_mag << std::endl;
        }
#endif
        if (tmp_hi_mag>m_hi_mag_double.local()) m_hi_mag_double = tmp_hi_mag;
    }
}

void MagnitudeLogger::synchronize(size_t icycle) {
    if (!m_input.static_tau) {
        m_enough_singles_for_dynamic_tau.update(icycle, m_nsingle>m_input.nenough_spawns_for_dynamic_tau);
        m_enough_doubles_for_dynamic_tau.update(icycle, m_nsingle>m_input.nenough_spawns_for_dynamic_tau);
        if (m_enough_singles_for_dynamic_tau && m_enough_doubles_for_dynamic_tau) {
            m_hi_mag_single.mpi_max();
            m_hi_mag_double.mpi_max();
            auto hi_mag_sum = m_hi_mag_single.reduced() + m_hi_mag_double.reduced();
            ASSERT(hi_mag_sum > 0.0);
            /*
             * highest transferred weight ~ tau x max(helem/prob)
             * i.e. recommended tau = max bloom / max(helem/prob)
             */
            m_tau = m_input.max_bloom / hi_mag_sum;
            // go halfway to predicted value
            m_psingle = (m_psingle+m_hi_mag_single / hi_mag_sum)/2;
            m_psingle = std::max(m_psingle, m_input.min_excit_class_prob);
            m_hi_mag_single = std::numeric_limits<defs::ham_comp_t>::min();
            m_hi_mag_double = std::numeric_limits<defs::ham_comp_t>::min();
        }
    }
}

MagnitudeLogger::MagnitudeLogger(const Options &input, size_t nsite, size_t nelec) :
        MagnitudeLogger(input, input.psingle_initial ? input.psingle_initial : psingle_guess(nsite, nelec)){}
