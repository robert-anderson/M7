//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/NumericField.h>
#include "gtest/gtest.h"
#include "src/core/list/List.h"

struct TestList : public List{
    NumericField<int> counter;
    TestList(size_t ncount=1, size_t nsegment=1):
    List("test list", nsegment),
    counter(this, ncount)
    {}
};


TEST(List, ThreadSafety) {
    size_t nrow = 3600;
    TestList list;
    list.expand(nrow);

#pragma omp parallel for default(none), shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push();
        list.counter(irow) = irow;
    }

    for (auto irow{0ul}; irow < list.nrow_per_segment(); ++irow) {
        ASSERT_EQ(list.counter(irow), irow);
    }
}

TEST(List, ThreadSerialization) {
    /*
     * Pushing to a List concurrently involves the resolution of thread contention
     * since the pointer to the next available position in the list must be atomically
     * incremented. This test is to make sure this contention is being resolved without
     * effective loss of thread parallelism
     */
    const size_t nrow = 1e7;
    TestList list;
    list.expand(nrow);

    auto do_something = [](size_t i){while (i%213){i = 3*i+1;}};

#pragma omp parallel for default(none), shared(list, do_something)
    for (size_t i = 0; i < nrow; ++i) {
        do_something(i);
        size_t irow = list.push();
        list.counter(irow) = irow;
    }
}

TEST(List, Communication){
    const size_t nrow = 100;
    const size_t nint = 6;
    TestList recv(nint);
    TestList send(nint, mpi::nrank());
    send.recv(&recv);
    send.expand(nrow);

    auto flat_index = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                         const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return (((isend * nrow + irecv) * nrow + irow * nrow) + ientry) % modular_divisor;
    };

    for (auto idst_rank{0ul}; idst_rank < mpi::nrank(); ++idst_rank) {
        for (auto irow{0ul}; irow < nrow; ++irow) {
            ASSERT_EQ(send.push(idst_rank), irow);
            for (auto iint{0ul}; iint < nint; ++iint) {
                send.counter(irow, idst_rank, iint) = flat_index(mpi::irank(), idst_rank, irow, iint);
            }
        }
    }

    send.communicate();
    /*
     * check that all entries came through as expected
     */
    size_t
    irow_tot = 0;
    for (auto isend{0ul}; isend < mpi::nrank(); ++isend) {
        for (auto irow{0ul}; irow < nrow; ++irow) {
            for (auto iint{0ul}; iint < nint; ++iint) {
                ASSERT_EQ(
                    recv.counter(irow_tot, 0, iint),
                    flat_index(isend, mpi::irank(), irow, iint)
                );
            }
            irow_tot++;
        }
    }
}

TEST(List, AllGather) {
    const size_t nrow = 10;
    const size_t nint = 1;
    TestList part(nint);
    TestList full(nint);
    for (size_t irow=0ul; irow<nrow; ++irow){
        ASSERT_EQ(part.expand_push(), irow);
        part.counter(irow)=irow+nrow*mpi::irank();
    }
    full.all_gather(part);
    full.print();
    for (size_t irow=0ul; irow<nrow*mpi::nrank(); ++irow){
        ASSERT_EQ(*full.counter(irow), irow);
    }
}

TEST(List, AllGatherRagged) {
    const size_t salt = 213;
    const size_t mod = 10;
    const size_t nrow = (salt*(mpi::irank()+1))%mod;
    const size_t nint = 1;
    TestList part(nint);
    TestList full(nint);
    for (size_t irow=0ul; irow<nrow; ++irow){
        ASSERT_EQ(part.expand_push(), irow);
        part.counter(irow)=mpi::irank()*mod+irow;
    }
    full.all_gather(part);
    size_t iflat = 0ul;
    for (size_t isrc=0ul; isrc<mpi::nrank(); ++isrc){
        for (size_t irow=0ul; irow<(salt*(isrc+1))%mod; ++irow){
            ASSERT_EQ(*full.counter(iflat), isrc*mod+irow);
            ++iflat;
        }
    }
}