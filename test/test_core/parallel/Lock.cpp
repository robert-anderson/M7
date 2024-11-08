//
// Created by Robert John Anderson on 06/09/2023.
//

#include "test_core/defs.h"
#include "M7_lib/parallel/SharedArray.h"
#include "M7_lib/parallel/Lock.h"

TEST(Lock, Basic) {
    lock::Vector lock_vec;
    const uint_t nelement = 1ul;
    lock_vec.resize(nelement);
    SharedArray<uint_t> array(nelement);
    const uint_t nterm=10ul;
    for (uint_t iwrite = 0ul; iwrite < nterm; ++iwrite) {
        lock_vec.acquire(0);
        const auto current = array[0];
        array.set_(0, current + iwrite);
        lock_vec.free(0);
    }
    mpi::barrier();
    const auto one_rank_tot = (nterm*(nterm-1))/2;
    ASSERT_EQ(array[0], mpi::nrank() * one_rank_tot);
}