//
// Created by rja on 15/08/2021.
//

#ifndef M7_MAES_H
#define M7_MAES_H

#include "src/core/bilinear/Bilinears.h"
#include "src/core/observables/RefExcits.h"

struct Maes {
    Epoch m_accum_epoch;
    Bilinears m_bilinears;
    RefExcits m_ref_excits;
    const size_t m_period;
    size_t m_icycle_period_start = ~0ul;

    Maes(const fciqmc_config::AvEsts &opts, const Hamiltonian &ham) :
            m_accum_epoch("MAE accumulation"), m_bilinears(opts, ham, m_accum_epoch),
            m_ref_excits(opts.m_ref_excits, ham.nsite()), m_period(opts.m_stats_period) {}

    operator bool() const {
        return m_bilinears || m_ref_excits;
    }

    bool all_stores_empty() const {
        return m_bilinears.all_stores_empty() && m_ref_excits.all_stores_empty();
    }

    bool is_period_cycle(size_t icycle) {
        if (!m_accum_epoch) return false;
        if (!m_period) return false;
        if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
            m_icycle_period_start = icycle;
            return false;
        }
        return !((icycle - m_icycle_period_start) % m_period);
    }

    void end_cycle() {
        m_bilinears.end_cycle();
    }
};


#endif //M7_MAES_H
