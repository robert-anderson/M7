//
// Created by rja on 12/12/22.
//

#include "InstEsts.h"

InstEsts::InstEsts(const sys::Sector sector, const References* refs, const conf::InstEsts& opts) {
    REQUIRE_TRUE(refs, "Invalid pointer to reference MBFs");
    if (opts.m_spin_square.m_value)
        m_spin_square = ptr::smart::make_unique<commuting_obs::SpinSquare>(sector.m_frm, refs);
}

void InstEsts::begin_cycle(uint_t icycle) {
    if (m_spin_square) m_spin_square->m_est.begin_cycle(icycle);
}

void InstEsts::make_numerator_contribs(const Walker& walker) {
    if (m_spin_square) m_spin_square->m_est.make_numerator_contribs(walker);
}

void InstEsts::end_cycle(uint_t icycle) {
    if (m_spin_square) m_spin_square->m_est.end_cycle(icycle);
}