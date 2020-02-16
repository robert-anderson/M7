//
// Created by Robert John Anderson on 2020-02-14.
//

#include <gtest/gtest.h>
#include <src/parallel/MPIWrapper.h>
#include <src/data/Specification.h>
#include <src/enumerators/CombinationEnumerator.h>
#include "src/data/PerforableMappedTable.h"


TEST(PerforableMappedTable, DataIntegrity) {
    const size_t nspatorb = 6;
    const size_t nelec = 4;

    Specification specification;
    specification.add<Determinant>(nspatorb);
    specification.add<size_t>(1);
    const size_t nrow = integer_utils::combinatorial(2*nspatorb, nelec);
    PerforableMappedTable<Determinant> table(specification, nrow);

    Determinant det(nspatorb);
    CombinationEnumerator enumerator(nspatorb*2, nelec);
    defs::inds inds(nelec);
    size_t count{~0ul};
    while(enumerator.next(inds, count)){
        det.set(inds);
        auto irow = table.push(0, det);
        *table.view<size_t>(irow) = count;
    }

    size_t irow;
    det.set(defs::inds{0, 6, 8, 11});
    irow = table.remove(0, det);
    ASSERT_EQ(irow, 151);
    ASSERT_THROW(table.lookup_view<Determinant>(det), KeyException);

    ASSERT_TRUE(table.lookup_view<Determinant>(det).is_zero());
    det.set(defs::inds{0, 1, 2, 3, 4});
    irow = table.push(0, det);
    ASSERT_EQ(irow, 151);
}