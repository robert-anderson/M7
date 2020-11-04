//
// Created by rja on 02/07/2020.
//

#include "ParallelizationStatsFile.h"

#if 0
ParallelizationStatsFile::ParallelizationStatsFile(const Options &input) :
    StatsFile(input.parallel_stats_path + '.' + std::to_string(mpi::irank())),
    m_cycle_number(this, 1, "Cycle number"),
    m_synchronization_wait_time(this, 1, "Seconds of idle time due to process sychronization"),
    m_nrow_free_walker_list(this, 1, "Number of free rows in walker list below high water mark"),
    m_walker_list_high_water_mark(this, 1, "High water mark of walker list"),
    m_walker_list_high_water_mark_fraction(this, 1, "High water mark of walker list as a fraction of current capacity"),
    m_nrow_sent(this, 1, "Total number of rows sent"),
    m_largest_nrow_sent(this, 1, "Largest number of rows sent to another rank"),
    m_largest_send_list_filled_fraction(this, 1,
                                        "Largest number of rows sent to another rank as a fraction of the segment capacity"),
    m_irank_largest_nrow_sent(this, 1, "MPI_COMM_WORLD index of the rank to which the largest number of rows was sent"),
    m_nrow_recv(this, 1, "Total number of rows received"),
    m_recv_list_filled_fraction(this, 1, "Total number of rows received as a fraction of current capacity") {}

#endif