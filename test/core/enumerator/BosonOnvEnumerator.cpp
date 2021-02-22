//
// Created by rja on 05/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/enumerator/BosonOnvEnumerator.h"
#include "src/core/field/BufferedFields.h"

TEST(BosonOnvEnumerator, Test){
    const size_t nmode=6ul, occ_cutoff=4;
    BosonOnvEnumerator enumerator(nmode, occ_cutoff);
    size_t i=~0ul;
    buffered::BosonOnv bonv(nmode);
    while (enumerator.next(bonv, i)){
        std::cout << bonv.to_string() << std::endl;
    }
    ASSERT_EQ(i, ci_utils::boson_dim(nmode, occ_cutoff));
}