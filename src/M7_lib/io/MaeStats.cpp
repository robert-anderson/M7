//
// Created by Robert J. Anderson on 18/08/2021.
//

#include "MaeStats.h"

MaeStatsRow::MaeStatsRow(bool rdms) :
        m_icycle(this, "Cycle number", false),
        m_total_norm((rdms ? this : nullptr), "Total WF norm estimate", true, {true, 9}),
        m_rdm_energy((rdms ? this : nullptr), "Energy estimate from RDMs", true, {true, 9}){
    DEBUG_ASSERT_TRUE(m_size, "row should have non-zero size");
}
