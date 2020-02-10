//
// Created by Robert John Anderson on 2020-02-09.
//

#include <gtest/gtest.h>
#include <src/parallel/MPIWrapper.h>
#include "src/data/BitfieldHasher.h"
#include "src/data/Specification.h"
#include "src/data/DataTable.h"

TEST(DataTable, AllToAllV) {

    //Specification spec;
    std::cout << numtypes::itype<float> << std::endl;

    exit(0);
    //MPIWrapper mpi;
    std::array<size_t, nnumeric> numeric_lengths;
    numeric_lengths[type_number<float>] = 2;
    defs::inds bitfield_lengths{34, 12};
    //DataTable send(numeric_lengths, bitfield_lengths, 4, mpi.nrank());
    DataTable send(numeric_lengths, bitfield_lengths, 4, 1);
    //send.print();
}