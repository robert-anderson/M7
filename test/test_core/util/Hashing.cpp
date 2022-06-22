//
// Created by Robert John Anderson on 2020-04-01.
//

#include <M7_lib/util/Hash.h>
#include "gtest/gtest.h"

TEST(Hashing, Determinism){
    const size_t nattempt = 10;

    hash::digest_t benchmark;
    size_t key = 99ul;
    benchmark = hash::fnv((char*)&key, sizeof(key));
    for (size_t i=0ul; i<nattempt; ++i){
        ASSERT_EQ(hash::fnv((char*)&key, sizeof(key)), benchmark);
    }
}