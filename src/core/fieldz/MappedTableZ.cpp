//
// Created by rja on 18/02/2021.
//

#include "MappedTableZ.h"

MappedTableBaseZ::MappedTableBaseZ(size_t nbucket) : m_buckets(nbucket){}

size_t MappedTableBaseZ::nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket) {
    // number of skips is on average half the average number of rows per bin.
    return std::max(c_nbucket_min, size_t(nrow_max / double(2 * mean_skips_per_bucket)));
}

size_t MappedTableBaseZ::nbucket() const {
    return m_buckets.size();
}

constexpr double MappedTableBaseZ::c_default_remap_ratio;
constexpr size_t MappedTableBaseZ::c_default_remap_nlookup;
constexpr size_t MappedTableBaseZ::c_nbucket_min;
