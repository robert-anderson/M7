//
// Created by Robert John Anderson on 2020-02-23.
//

#include <gtest/gtest.h>
#include <src/defs.h>
#include <src/io/FileIterator.h>
#include "src/io/StatsFile.h"

TEST(StatsFile, CorrectFormat){
    StatsFile stats_file("M7.stats");
    stats_file.add_column<size_t>("Iteration number");
    stats_file.add_column<double>("Shift");
    stats_file.add_column<size_t>("Number of Occupied Determinants");
    stats_file.add_column<defs::ham_t>("Reference Projected Energy Numerator");
    stats_file.add_column<defs::ham_t>("Reference Projected Energy Denominator");
    stats_file.add_column<size_t>("Number of Initiators");
    stats_file.write_header();

    FileIterator file_iterator("M7.stats");
    std::string line;
    while (file_iterator.next(line)){
        std::cout << line << std::endl;
    }
}