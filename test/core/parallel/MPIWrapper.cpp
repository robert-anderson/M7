//
// Created by Robert John Anderson on 2020-04-15.
//

#include "gtest/gtest.h"
#include "src/core/parallel/MPIWrapper.h"


TEST(MPIWrapper, Alltoall){
    defs::inds send(mpi::nrank(), 0ul);
    defs::inds recv(mpi::nrank(), 0ul);
    const size_t salt = 213;
    for (size_t irecv=0ul; irecv<mpi::nrank(); ++irecv){
        send[irecv] = mpi::irank()+salt*(1+irecv);
    }
    mpi::all_to_all(send, recv);
    for (size_t isent=0ul; isent<mpi::nrank(); ++isent){
        ASSERT_EQ(recv[isent], isent+salt*(1+mpi::irank()));
    }
}