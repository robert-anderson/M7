//
// Created by rja on 13/04/2021.
//

#ifndef M7_PARALLELSTATS_H
#define M7_PARALLELSTATS_H

#include "src/core/field/Row.h"
#include "StatsTable.h"

struct ParallelStatsRow : Row {
    fields::Number<size_t> m_icycle;
    fields::Number<double> m_synchronization_overhead;
    fields::Number<size_t> m_nblock_wf_ra;
    fields::Number<defs::wf_comp_t> m_nwalker_total;
    fields::Number<size_t> m_nwalker_lookup_skip;
    fields::Number<size_t> m_nwalker_lookup;
    fields::Number<size_t> m_nrow_free_walker_list;
    fields::Number<size_t> m_walker_list_high_water_mark;
    fields::Number<double> m_walker_list_high_water_mark_fraction;
    fields::Number<size_t> m_nrow_sent;
    fields::Number<size_t> m_largest_nrow_sent;
    fields::Number<double> m_largest_send_list_filled_fraction;
    fields::Number<size_t> m_irank_largest_nrow_sent;
    fields::Number<size_t> m_nrow_recv;
    fields::Number<double> m_recv_list_filled_fraction;

    ParallelStatsRow();
};

typedef StatsTable<ParallelStatsRow> ParallelStats;

#endif //M7_PARALLELSTATS_H
