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
    StatsField<defs::ham_comp_t> m_nwalker;
    StatsField<size_t> m_noccupied_det;
    StatsField<double> m_synchronization_wait_time;

    explicit ParallelizationStatsFile(const Options &input);
};



#endif //M7_PARALLELIZATIONSTATSFILE_H
