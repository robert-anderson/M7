//
// Created by rja on 02/07/2020.
//

#include "ParallelizationStatsFile.h"

ParallelizationStatsFile::ParallelizationStatsFile(const Options &input) :
        StatsFile(input.parallel_stats_path+'.'+std::to_string(mpi::irank())),
        m_cycle_number(this, 1, "Cycle number"),
        m_nwalker(this, 1, "Total number of walkers"),
        m_noccupied_det(this, 1, "Number of occupied determinants"),
        m_synchronization_wait_time(this, 1, "Seconds of idle time due to process sychnonization")
{}
