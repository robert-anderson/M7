//
// Created by rja on 18/08/2021.
//

#ifndef M7_MAESTATS_H
#define M7_MAESTATS_H

#include "field/Row.h"
#include "StatsTable.h"

struct MaeStatsRow : StatsRow {
    statistic::Number<size_t> m_icycle;
    statistic::Number<defs::wf_t> m_total_norm;
    statistic::Number<defs::ham_comp_t> m_rdm_energy;

    MaeStatsRow(bool rdms, bool spec_moms);
};

typedef StatsTable<MaeStatsRow> MaeStats;


#endif //M7_MAESTATS_H
