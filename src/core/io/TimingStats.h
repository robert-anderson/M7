//
// Created by rja on 23/07/2021.
//

#ifndef M7_TIMINGSTATS_H
#define M7_TIMINGSTATS_H

#include "field/Row.h"
#include "StatsTable.h"

struct TimingStatsRow : StatsRow {
    statistic::Number<double> m_total_synchronization_overhead;
    statistic::Number<double> m_propagate_loop_time;
    statistic::Number<double> m_communication_time;
    statistic::Number<double> m_annihilation_loop_time;
    statistic::Number<double> m_total_cycle_time;

    TimingStatsRow() :
        m_total_synchronization_overhead(this, "Total time waited at MPI_Barrier"),
        m_propagate_loop_time(this, "Time spent in loop over occupied ONVs"),
        m_communication_time(this, "Time spent in communicating spawns"),
        m_annihilation_loop_time(this, "Time spent in annihilation loop"),
        m_total_cycle_time(this, "Total cycle time"){}
};

typedef StatsTable<TimingStatsRow> TimingStats;

#endif //M7_TIMINGSTATS_H
