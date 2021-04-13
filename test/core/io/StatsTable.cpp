//
// Created by rja on 12/04/2021.
//

#include "gtest/gtest.h"
#include "src/core/io/FciqmcStats.h"

TEST(StatsTable, Test) {
    FciqmcStats stats_table("test.stats", "FCIQMC", {1, 1});
    stats_table.m_row.m_tau = 123;
    stats_table.flush();
}
