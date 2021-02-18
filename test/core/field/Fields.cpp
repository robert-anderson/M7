//
// Created by rja on 26/01/2021.
//

#include <src/core/sample/PRNG.h>
#include <src/core/fieldz/BufferedFields.h>
#include "gtest/gtest.h"

TEST(Fields, HashUniformityTrueRandom){
    const size_t nsite = 50;
    const size_t nelec = 28;
    const size_t nbucket = 20;
    const size_t ndraw = 10000000;

    defs::inds freqs(nbucket, 0ul);
    PRNG prng(14, 10000);

    buffered::FermionOnv fonv(nsite);

    for (size_t idraw=0ul; idraw<ndraw; ++idraw){
        fonv.zero();
        for (auto ispin: {0ul, 1ul}) {
            size_t nset = 0ul;
            while (nset < nelec/2) {
                size_t isite = prng.draw_uint(nsite);
                if (fonv.get(ispin, isite)) continue;
                fonv.set(ispin, isite);
                ++nset;
            }
        }
        auto ibucket = fonv.hash()%nbucket;
        freqs[ibucket]++;
    }
    auto tot = std::accumulate(freqs.begin(), freqs.end(), 0);
    // check that the hash is uniform within 1%
    for (const auto& freq: freqs) ASSERT_LT(std::abs(1.0-nbucket*freq/double(tot)), 0.01);
}


TEST(Fields, HashUniformityLowIndexMoreLikely){
    /*
     * simulates the occupation of a CI vector
     */
    const size_t nsite = 50;
    const size_t nelec = 28;
    const size_t nbucket = 20;
    const size_t ndraw = 10000000;
    const defs::prob_t prob = 0.6;

    defs::inds freqs(nbucket, 0ul);
    PRNG prng(14, 10000);

    buffered::FermionOnv fonv(nsite);

    for (size_t idraw=0ul; idraw<ndraw; ++idraw){
        fonv.zero();
        for (auto ispin: {0, 1}) {
            size_t nset = 0ul;
            size_t isite = 0ul;
            while (nset<nelec/2){
                if (!fonv.get(ispin, isite)) {
                    if (prng.draw_float() < prob) {
                        fonv.set(ispin, isite);
                        ++nset;
                    }
                }
                ++isite;
                isite%=nsite;
            }
        }
        auto ibucket = fonv.hash()%nbucket;
        freqs[ibucket]++;
    }
    auto tot = std::accumulate(freqs.begin(), freqs.end(), 0);
    // check that the hash is uniform within 1%
    for (const auto& freq: freqs) ASSERT_LT(std::abs(1.0-nbucket*freq/double(tot)), 0.01);
}