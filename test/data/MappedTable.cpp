//
// Created by Robert John Anderson on 2020-02-13.
//

#include <gtest/gtest.h>
#include <src/parallel/MPIWrapper.h>
#include <src/data/Specification.h>
#include <src/enumerators/CombinationEnumerator.h>
#include "src/data/MappedTable.h"

/*
TEST(MappedTable, DataIntegrity) {
    const size_t nspatorb = 6;
    const size_t nelec = 4;

    Specification specification;
    specification.add<Determinant>(nspatorb);
    specification.add<size_t>(1);
    const size_t nrow = integer_utils::combinatorial(2*nspatorb, nelec);
    MappedTable<Determinant> table(specification, nrow);

    Determinant det(nspatorb);
    CombinationEnumerator enumerator(nspatorb*2, nelec);
    defs::inds inds(nelec);
    size_t count{~0ul};
    while(enumerator.next(inds, count)){
        det.set(inds);
        auto irow = table.push(0, det);
        *table.view<size_t>(irow) = count;
    }
    ASSERT_TRUE(table.highwatermark()[0]==nrow);
    det.set(defs::inds{0, 6, 8, 11});
    ASSERT_EQ(table.lookup(det, 0), 151);
    det.set(defs::inds{8, 9, 10, 11});
    ASSERT_EQ(table.lookup(det, 0), nrow-1);
    det.set(defs::inds{7, 8, 9, 10, 11});
    ASSERT_EQ(table.lookup(det, 0), ~0ul);
}


TEST(MappedTable, ThreadSafety) {
    const size_t nspatorb = 6;
    const size_t nelec = 4;

    Specification specification;
    specification.add<Determinant>(nspatorb);
    specification.add<size_t>(1);
    const size_t nrow = integer_utils::combinatorial(2 * nspatorb, nelec);
    MappedTable<Determinant> table(specification, nrow);

    auto all_combs = CombinationEnumerator(nspatorb * 2, nelec).enumerate();
#pragma omp parallel
    {
        Determinant det(nspatorb);
#pragma omp for
        for (size_t icomb = 0ul; icomb < nrow; ++icomb) {
            det.set(all_combs[icomb]);
            auto mutex = table.get_mutex()
            auto irow = table.push(0, det);
            *table.view<size_t>(irow) = icomb;
        }
    }

    Determinant det(nspatorb);
    ASSERT_EQ(table.highwatermark()[0], nrow);

    det.set(defs::inds{7, 8, 9, 10, 11});
    ASSERT_EQ(table.lookup(det, 0), ~0ul);
}
 */