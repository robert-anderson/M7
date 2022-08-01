//
// Created by Robert John Anderson on 2020-04-01.
//

#include <M7_lib/util/Hash.h>
#include "gtest/gtest.h"

TEST(Hashing, Determinism){
    const uint_t nattempt = 10;

    hash::digest_t benchmark;
    uint_t key = 99ul;
    benchmark = hash::fnv((buf_t *)&key, sizeof(key));
    for (uint_t i=0ul; i<nattempt; ++i){
        ASSERT_EQ(hash::fnv((buf_t *)&key, sizeof(key)), benchmark);
    }
}