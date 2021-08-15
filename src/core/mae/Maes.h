//
// Created by rja on 15/08/2021.
//

#ifndef M7_MAES_H
#define M7_MAES_H

#include "src/core/bilinear/BilinearEstimators.h"
#include "src/core/observables/RefExcits.h"

struct Maes {
    BilinearEstimators m_bilinears;
    RefExcits m_ref_excits;
    const size_t m_period;
    size_t m_icycle_period_start = ~0ul;

    Maes(const fciqmc_config::AvEsts& opts, const Hamiltonian& ham):
        m_bilinears(opts, ham), m_ref_excits(opts.m_ref_excits, ham.nsite()),
        m_period(opts.m_stats_period);

    bool is_period_cycle(size_t icycle) {
        if (!m_accum_epoch) return false;
        if (!m_period) return false;
        if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
            m_icycle_period_start = icycle;
            return false;
        }
        return !((icycle - m_icycle_period_start) % m_period);
    }

};


#endif //M7_MAES_H
