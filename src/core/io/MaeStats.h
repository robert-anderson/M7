//
// Created by rja on 18/08/2021.
//

#ifndef M7_MAESTATS_H
#define M7_MAESTATS_H

#include "src/core/field/Row.h"
#include "StatsTable.h"

struct MaeStatsRow : Row {
    field::Number<size_t> m_icycle;
    field::Number<defs::wf_t> m_total_norm;
    field::Number<defs::ham_comp_t> m_rdm_energy;

    MaeStatsRow(bool rdms, bool spec_moms, bool have_energy);
};

typedef StatsTable<MaeStatsRow> MaeStats;


#endif //M7_MAESTATS_H
