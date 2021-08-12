//
// Created by rja on 13/04/2021.
//

#ifndef M7_PARALLELSTATS_H
#define M7_PARALLELSTATS_H

#include "src/core/field/Row.h"
#include "StatsTable.h"

struct ParallelStatsRow : Row {
    field::Number<size_t> m_icycle;
    field::Number<double> m_synchronization_overhead;
    field::Number<size_t> m_nblock_wf_ra;
    field::Number<defs::wf_comp_t> m_nwalker_total;
    field::Number<size_t> m_nwalker_lookup_skip;
    field::Number<size_t> m_nwalker_lookup;
    field::Number<size_t> m_nrow_free_walker_list;
    field::Number<size_t> m_walker_list_high_water_mark;
    field::Number<double> m_walker_list_high_water_mark_fraction;
    field::Number<size_t> m_nrow_sent;
    field::Number<size_t> m_largest_nrow_sent;
    field::Number<double> m_largest_send_list_filled_fraction;
    field::Number<size_t> m_irank_largest_nrow_sent;
    field::Number<size_t> m_nrow_recv;
    field::Number<double> m_recv_list_filled_fraction;

    ParallelStatsRow();
};

typedef StatsTable<ParallelStatsRow> ParallelStats;

#endif //M7_PARALLELSTATS_H
