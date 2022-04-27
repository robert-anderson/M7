//
// Created by Robert J. Anderson on 18/02/2021.
//

#include "MappedTable.h"

MappedTableBase::MappedTableBase(size_t nbucket, size_t remap_nlookup, double remap_ratio) :
    m_buckets(nbucket?nbucket:c_default_nbucket),
    m_remap_nlookup(remap_nlookup?remap_nlookup:c_default_remap_nlookup),
    m_remap_ratio(remap_ratio!=0.0?remap_ratio:c_default_remap_ratio){}

size_t MappedTableBase::nbucket_guess(size_t nitem, double ratio) {
    return nitem / (2*ratio + 1);
}

size_t MappedTableBase::nbucket() const {
    return m_buckets.size();
}

void MappedTableBase::clear_map() {
    for (auto &bucket: m_buckets) bucket.clear();
}

bool MappedTableBase::remap_due() const {
    return (m_nlookup_total >= m_remap_nlookup) &&
           (double(m_nskip_total) / double(m_nlookup_total)) > m_remap_ratio;
}

bool MappedTableBase::all_nonzero_rows_mapped(const TableBase &source) const {
    // make a set out of all bucketed indices
    std::set<size_t> set;
    for (const auto& bucket: m_buckets){
        for (const auto& item: bucket) set.insert(item);
    }
    // then go through table, if the row is non-zero, its index should be in the set.
    for (size_t irow=0ul; irow<source.m_hwm; ++irow){
        bool in_set = set.find(irow)!=set.end();
        if (source.is_cleared(irow)==in_set) return false;
    }
    return true;
}

constexpr size_t MappedTableBase::c_default_nbucket;
constexpr double MappedTableBase::c_default_remap_ratio;
constexpr size_t MappedTableBase::c_default_remap_nlookup;