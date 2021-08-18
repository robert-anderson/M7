//
// Created by rja on 18/08/2021.
//

#include "MaeStats.h"

MaeStatsRow::MaeStatsRow(bool rdms, bool spec_moms, bool have_energy) :
        m_icycle(this, "Cycle number"),
        m_total_norm((rdms || spec_moms ? this : nullptr), "Total WF norm estimate"),
        m_rdm_energy((rdms ? this : nullptr), "Energy estimate from RDMs"){}
