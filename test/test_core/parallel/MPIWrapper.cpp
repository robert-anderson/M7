//
// Created by Robert John Anderson on 2020-04-15.
//

#include <M7_lib/util/Hash.h>
#include <algorithm>
#include "gtest/gtest.h"
#include "M7_lib/parallel/MPIWrapper.h"

TEST(MPIWrapper, AllSum){
    ASSERT_EQ(mpi::nrank(), mpi::all_sum(1ul));
    uint_t i = mpi::irank()+1;
    uint_t res = mpi::all_sum(i);
    ASSERT_EQ(res, (mpi::nrank()*(mpi::nrank()+1))/2);
}

TEST(MPIWrapper, AllMax){
    uint_t i = hash::in_range(mpi::irank(), 5, 19);
    uint_t res = mpi::all_max(i);
    uintv_t chk(mpi::nrank());
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) chk[irank] = hash::in_range(irank, 5, 19);
    std::sort(chk.begin(), chk.end());
    ASSERT_EQ(res, chk.back());
}

TEST(MPIWrapper, AllMin){
    uint_t i = hash::in_range(mpi::irank(), 5, 19);
    uint_t res = mpi::all_min(i);
    uintv_t chk(mpi::nrank());
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) chk[irank] = hash::in_range(irank, 5, 19);
    std::sort(chk.begin(), chk.end());
    ASSERT_EQ(res, chk.front());
}

TEST(MPIWrapper, Alltoall){
    uintv_t send(mpi::nrank(), 0ul);
    uintv_t recv(mpi::nrank(), 0ul);
    for (uint_t irecv=0ul; irecv<mpi::nrank(); ++irecv){
        send[irecv] = hash::in_range({irecv, mpi::irank()}, 3, 123);
    }
    mpi::all_to_all(send, recv);
    for (uint_t isent=0ul; isent<mpi::nrank(); ++isent){
        ASSERT_EQ(recv[isent], hash::in_range({mpi::irank(), isent}, 3, 123));
    }
}

TEST(MPIWrapper, Allgatherv){
    const uint_t n=4;
    uintv_t send(n, 0ul);
    uintv_t recv(n * mpi::nrank(), 0ul);
    for (uint_t i=0ul; i<n; ++i){
        send[i] = hash::in_range({mpi::irank(), i}, 3, 123);
    }
    uintv_t recvcounts(mpi::nrank(), n);
    uintv_t displs(mpi::nrank(), 0);
    for (uint_t i=1ul; i<mpi::nrank(); ++i) displs[i] = displs[i-1]+n;

    mpi::all_gatherv(send.data(), n, recv.data(), recvcounts, displs);
    uint_t iflat = 0ul;
    for (uint_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (uint_t i=0ul; i<n; ++i){
            ASSERT_EQ(recv[iflat], hash::in_range({isrc, i}, 3, 123));
            ++iflat;
        }
    }
}

TEST(MPIWrapper, Somegatherv){
    const uint_t n=4;
    uintv_t send(n, 0ul);
    uintv_t recv(n * mpi::nrank(), 0ul);
    for (uint_t i=0ul; i<n; ++i){
        send[i] = hash::in_range({mpi::irank(), i}, 3, 123);
    }
    // suppose that only even-indexed ranks gather data
    auto recving = [](uint_t irank=mpi::irank()) -> bool {return !(irank&1ul);};
    uintv_t sendcounts(mpi::nrank(), 0ul);
    uintv_t senddispls(mpi::nrank(), 0ul);
    uintv_t recvcounts(mpi::nrank(), 0ul);
    uintv_t recvdispls(mpi::nrank(), 0ul);
    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank){
        // all ranks send all their data to only the receiving ranks
        if (recving(irank)) sendcounts[irank] = n;
    }
    if (recving()) recvcounts.assign(mpi::nrank(), n);
    recvdispls = mpi::counts_to_displs_consec(recvcounts);

    mpi::all_to_allv(send.data(), sendcounts, senddispls, recv.data(), recvcounts, recvdispls);
    uint_t iflat = 0ul;
    for (uint_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (uint_t i=0ul; i<n; ++i){
            if (mpi::irank()%2) ASSERT_EQ(recv[iflat], 0);
            else ASSERT_EQ(recv[iflat], hash::in_range({isrc, i}, 3, 123));
            ++iflat;
        }
    }
}

TEST(MPIWrapper, AllgathervRagged){
    const uint_t nsend=hash::in_range(mpi::irank(), 5, 17);
    uintv_t send(nsend, 0ul);
    for (uint_t i=0ul; i<nsend; ++i){
        send[i] = hash::in_range({mpi::irank(), i}, 3, 123);
    }

    uintv_t recvcounts(mpi::nrank(), 0);
    mpi::all_gather(nsend, recvcounts);
    for (uint_t i=0ul; i<mpi::nrank(); ++i){
        ASSERT_EQ(recvcounts[i], hash::in_range(i, 5, 17));
    }

    uintv_t displs(mpi::nrank(), 0);
    for (uint_t i=1ul; i<mpi::nrank(); ++i) displs[i] = displs[i-1]+recvcounts[i-1];

    const uint_t nrecv = displs.back()+recvcounts.back();
    uintv_t recv(nrecv, 0ul);

    mpi::all_gatherv(send.data(), nsend, recv.data(), recvcounts, displs);

    uint_t iflat = 0ul;
    for (uint_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (uint_t i=0ul; i<recvcounts[isrc]; ++i){
            ASSERT_EQ(recv[iflat], hash::in_range({isrc, i}, 3, 123));
            ++iflat;
        }
    }
}

TEST(MPIWrapper, MaxLocMinLoc){
    typedef double T;
    auto to_T = [](uint_t irank){
        return std::pow(1.13, (irank+7*3)%5+1);
    };
    T local = to_T(mpi::irank());
    std::pair<T, uint_t> max{std::numeric_limits<T>::min(), 0};
    std::pair<T, uint_t> min{std::numeric_limits<T>::max(), 0};
    for (uint_t i=0ul; i<mpi::nrank(); ++i) {
        auto d = to_T(i);
        if (d>max.first) max = {d, i};
        if (d<min.first) min = {d, i};
    }

    std::pair<T, uint_t> mpi_max{std::numeric_limits<T>::min(), 0};
    std::pair<T, uint_t> mpi_min{std::numeric_limits<T>::max(), 0};

    ASSERT_EQ(MPI_DOUBLE_INT, mpi::pair_type<T>());

    ASSERT_TRUE(mpi::all_maxloc(local, mpi_max));
    ASSERT_EQ(max, mpi_max);
    ASSERT_TRUE(mpi::all_minloc(local, mpi_min));
    ASSERT_EQ(min, mpi_min);
}

TEST(MPIWrapper, RanksWithNonzero){
    uintv_t iranks_chk = {};
    for (uint_t irank=0ul; irank < mpi::nrank(); ++irank)
        if (hash::in_range(irank+123, 0, 2)) iranks_chk.push_back(irank);
    const bool i_have_nonzero = hash::in_range(mpi::irank()+123, 0, 2);
    const auto iranks = mpi::filter(i_have_nonzero);
    ASSERT_EQ(iranks, iranks_chk);
}

