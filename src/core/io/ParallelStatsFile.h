//
// Created by rja on 02/07/2020.
//

#ifndef M7_PARALLELSTATSFILE_H
#define M7_PARALLELSTATSFILE_H


#include "StatsFile.h"

struct ParallelStatsSpecifier : StatsSpecifier {
    StatsColumn<size_t> m_icycle;
    StatsColumn<double> m_synchronization_overhead;
    StatsColumn<defs::wf_comp_t> m_nwalker;
    StatsColumn<size_t> m_nrow_free_walker_list;
    StatsColumn<size_t> m_walker_list_high_water_mark;
    StatsColumn<double> m_walker_list_high_water_mark_fraction;
    StatsColumn<size_t> m_nrow_sent;
    StatsColumn<size_t> m_largest_nrow_sent;
    StatsColumn<double> m_largest_send_list_filled_fraction;
    StatsColumn<size_t> m_irank_largest_nrow_sent;
    StatsColumn<size_t> m_nrow_recv;
    StatsColumn<double> m_recv_list_filled_fraction;

    ParallelStatsSpecifier() :
            StatsSpecifier("Parallelization"),
            m_icycle(this, "Cycle number"),
            m_synchronization_overhead(this, "Time waited at MPI_Barrier"),
            m_nwalker(this, "WF L1 norm (number of walkers)"),
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

#endif //M7_PARALLELSTATSFILE_H
