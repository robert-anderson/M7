//
// Created by rja on 02/07/2020.
//

#ifndef M7_PARALLELIZATIONSTATSFILE_H
#define M7_PARALLELIZATIONSTATSFILE_H

#include "StatsFile.h"
#include "Options.h"

struct ParallelizationStatsFile : public StatsFile {
public:
    StatsField<size_t> m_cycle_number;
    StatsField<double> m_synchronization_wait_time;
    StatsField<size_t> m_nrow_free_walker_list;
    StatsField<size_t> m_walker_list_high_water_mark;
    StatsField<double> m_walker_list_high_water_mark_fraction;
    StatsField<size_t> m_nrow_sent;
    StatsField<size_t> m_largest_nrow_sent;
    StatsField<double> m_largest_send_list_filled_fraction;
    StatsField<size_t> m_irank_largest_nrow_sent;
    StatsField<size_t> m_nrow_recv;
    StatsField<double> m_recv_list_filled_fraction;

    explicit ParallelizationStatsFile(const Options &input);
};



#endif //M7_PARALLELIZATIONSTATSFILE_H
