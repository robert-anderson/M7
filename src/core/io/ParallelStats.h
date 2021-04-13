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

    ParallelStatsRow() :
            m_icycle(this, "Cycle number"),
            m_synchronization_overhead(this, "Time waited at MPI_Barrier"),
            m_nblock_wf_ra(this, "Number of blocks of the WF rank allocator stored here"),
            m_nwalker_total(this, "Total WF L1 norm (number of walkers) over all WF parts"),
            m_nwalker_lookup_skip(this, "Number of skips required by walker hashtable lookups"),
            m_nwalker_lookup(this, "Number of walker hashtable lookups"),
            m_nrow_free_walker_list(this, "Free rows in the walker table"),
            m_walker_list_high_water_mark(this, "High water mark of walker table"),
            m_walker_list_high_water_mark_fraction(
                    this, "High water mark of water table as a fraction of allocated rows"),
            m_nrow_sent(this, "Total number of rows sent to other ranks"),
            m_largest_nrow_sent(this, "Largest number of rows sent to a single rank"),
            m_largest_send_list_filled_fraction(
                    this, "Largest number of rows sent to a single rank as a fraction of allocated rows"),
            m_irank_largest_nrow_sent(this, "Index of the rank to which the largest number of rows was sent"),
            m_nrow_recv(this, "Total number of rows received"),
            m_recv_list_filled_fraction(
                    this, "Total number of rows received as a fraction of allocated rows")
    {}
};

typedef StatsTable<ParallelStatsRow> ParallelStats;

#endif //M7_PARALLELSTATS_H
