//
// Created by Robert J. Anderson on 18/08/2021.
//

#ifndef M7_MAESTATS_H
#define M7_MAESTATS_H

#include <M7_lib/field/Row.h>

#include "StatsTable.h"

struct MaeStatsRow : StatsRow {
    statistic::Number<uint_t> m_icycle;
    statistic::Number<defs::wf_t> m_total_norm;
    statistic::Number<defs::ham_comp_t> m_rdm_energy;

    MaeStatsRow(bool rdms, bool spec_moms);
};

typedef StatsTable<MaeStatsRow> MaeStats;


#endif //M7_MAESTATS_H
