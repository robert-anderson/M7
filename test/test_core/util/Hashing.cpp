//
// Created by Robert John Anderson on 2020-04-01.
//

#include <M7_lib/util/Hashing.h>
#include "gtest/gtest.h"

TEST(Hashing, Determinism){
    const size_t nattempt = 10;

    utils::hash::digest_t benchmark;
    size_t key = 99ul;
    benchmark = utils::hash::fnv((char*)&key, sizeof(key));
    for (size_t i=0ul; i<nattempt; ++i){
        ASSERT_EQ(utils::hash::fnv((char*)&key, sizeof(key)), benchmark);
    }
}