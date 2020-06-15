//
// Created by Robert John Anderson on 2020-04-15.
//

#include "gtest/gtest.h"
#include "src/core/parallel/MPIWrapper.h"

TEST(MPIWrapper, AllSum){
    size_t i = mpi::irank()+1;
    size_t res = mpi::all_sum(i);
    ASSERT_EQ(res, (mpi::nrank()*(mpi::nrank()+1))/2);
}

TEST(MPIWrapper, AllMax){
    size_t i = (23*(mpi::irank()+1))%5;
    size_t res = mpi::all_max(i);
    size_t test = 0ul;
    for (size_t j=0ul; j<mpi::nrank(); ++j){
        auto tmp = (23*(j+1))%5;
        if (tmp>test) test = tmp;
    }
    ASSERT_EQ(res, test);
}

TEST(MPIWrapper, AllMin){
    size_t i = (23*(mpi::irank()+1))%5;
    size_t res = mpi::all_max(i);
    size_t test = 0ul;
    for (size_t j=0ul; j<mpi::nrank(); ++j){
        auto tmp = (23*(j+1))%5;
        if (tmp<test) test = tmp;
    }
    ASSERT_EQ(res, test);
}

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