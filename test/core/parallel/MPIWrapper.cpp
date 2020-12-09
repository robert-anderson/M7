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
    size_t res = mpi::all_min(i);
    size_t test = ~0ul;
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

TEST(MPIWrapper, Allgatherv){
    const size_t n=4;
    defs::inds send(n, 0ul);
    defs::inds recv(n*mpi::nrank(), 0ul);
    const size_t salt = 213;
    for (size_t i=0ul; i<n; ++i){
        send[i] = mpi::irank()+salt*(1+i);
    }
    defs::inds recvcounts(mpi::nrank(), n);
    defs::inds displs(mpi::nrank(), 0);
    for (size_t i=1ul; i<mpi::nrank(); ++i) displs[i] = displs[i-1]+n;

    mpi::all_gatherv(send.data(), n, recv.data(), recvcounts, displs);
    size_t iflat = 0ul;
    for (size_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (size_t i=0ul; i<n; ++i){
            ASSERT_EQ(recv[iflat], isrc+salt*(1+i));
            ++iflat;
        }
    }
}

TEST(MPIWrapper, AllgathervRagged){
    const size_t salt = 213;
    const size_t mod = 7;
    const size_t nsend=((mpi::irank()+1)*salt)%mod;
    defs::inds send(nsend, 0ul);
    for (size_t i=0ul; i<nsend; ++i){
        send[i] = mpi::irank()+salt*(1+i);
    }

    defs::inds recvcounts(mpi::nrank(), 0);
    mpi::all_gather(nsend, recvcounts);
    for (size_t i=0ul; i<mpi::nrank(); ++i){
        ASSERT_EQ(recvcounts[i], ((i+1)*salt)%mod);
    }

    defs::inds displs(mpi::nrank(), 0);
    for (size_t i=1ul; i<mpi::nrank(); ++i) displs[i] = displs[i-1]+recvcounts[i-1];

    const size_t nrecv = displs.back()+recvcounts.back();
    defs::inds recv(nrecv, 0ul);

    mpi::all_gatherv(send.data(), nsend, recv.data(), recvcounts, displs);

    size_t iflat = 0ul;
    for (size_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (size_t i=0ul; i<recvcounts[isrc]; ++i){
            ASSERT_EQ(recv[iflat], isrc+salt*(1+i));
            ++iflat;
        }
    }
}


TEST(MPIWrapper, MaxLocMinLoc){
    typedef double T;
    auto to_T = [](size_t irank){
        return std::pow(1.13, (irank+7*3)%5+1);
    };
    T local = to_T(mpi::irank());
    std::pair<T, size_t> max{std::numeric_limits<T>::min(), 0};
    std::pair<T, size_t> min{std::numeric_limits<T>::max(), 0};
    for (size_t i=0ul; i<mpi::nrank(); ++i) {
        auto d = to_T(i);
        if (d>max.first) max = {d, i};
        if (d<min.first) min = {d, i};
    }

    std::pair<T, size_t> mpi_max{std::numeric_limits<T>::min(), 0};
    std::pair<T, size_t> mpi_min{std::numeric_limits<T>::max(), 0};

    ASSERT_EQ(MPI_DOUBLE_INT, mpi_pair_type<T>());

    ASSERT_TRUE(mpi::all_maxloc(local, mpi_max));
    ASSERT_EQ(max, mpi_max);
    ASSERT_TRUE(mpi::all_minloc(local, mpi_min));
    ASSERT_EQ(min, mpi_min);
}

