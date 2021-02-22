//
// Created by rja on 18/02/2021.
//

#include "MappedTable.h"

MappedTableBase::MappedTableBase(size_t nbucket) : m_buckets(nbucket){}

size_t MappedTableBase::nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket) {
    // number of skips is on average half the average number of rows per bin.
    return std::max(c_nbucket_min, size_t(nrow_max / double(2 * mean_skips_per_bucket)));
}

size_t MappedTableBase::nbucket() const {
    return m_buckets.size();
}

constexpr double MappedTableBase::c_default_remap_ratio;
constexpr size_t MappedTableBase::c_default_remap_nlookup;
constexpr size_t MappedTableBase::c_nbucket_min;
