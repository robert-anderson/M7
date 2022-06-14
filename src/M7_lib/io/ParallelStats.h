//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_PARALLELSTATS_H
#define M7_PARALLELSTATS_H

#include <M7_lib/field/Row.h>

#include "StatsTable.h"

struct ParallelStatsRow : StatsRow {
    statistic::Number<size_t> m_icycle;
    statistic::Number<double> m_synchronization_overhead;
    statistic::Number<size_t> m_nblock_wf_ra;
    statistic::Number<defs::wf_comp_t> m_nwalker_total;
    statistic::Number<size_t> m_nwalker_lookup_skip;
    statistic::Number<size_t> m_nwalker_lookup;
    statistic::Number<size_t> m_nrow_free_walker_list;
    statistic::Number<size_t> m_walker_list_high_water_mark;
    statistic::Number<double> m_walker_list_high_water_mark_fraction;
    statistic::Number<size_t> m_nrow_sent;
    statistic::Number<size_t> m_largest_nrow_sent;
    statistic::Number<double> m_largest_send_list_filled_fraction;
    statistic::Number<size_t> m_irank_largest_nrow_sent;
    statistic::Number<size_t> m_nrow_recv;
    statistic::Number<double> m_recv_list_filled_fraction;

    ParallelStatsRow();
};

typedef StatsTable<ParallelStatsRow> ParallelStats;

#endif //M7_PARALLELSTATS_H
