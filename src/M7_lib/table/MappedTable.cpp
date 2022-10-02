//
// Created by Robert J. Anderson on 18/02/2021.
//

#include "MappedTable.h"

MappedTableBase::MappedTableBase(MappedTableOptions opts) : m_opts(opts) {}

MappedTableBase &MappedTableBase::operator=(const MappedTableBase &other) {
    if (this==&other) return *this;
    m_buckets = other.m_buckets;
    m_nskip_total = other.m_nskip_total;
    m_nlookup_total = other.m_nlookup_total;
    return *this;
}

MappedTableBase::MappedTableBase(const MappedTableBase &other) : MappedTableBase(other.m_opts){
    m_buckets = other.m_buckets;
}

bool MappedTableBase::operator==(const MappedTableBase &other) const {
    if (this==&other) return true;
    if (nbucket()!=other.nbucket()) return false;
    if (m_opts.m_remap_ratio!=other.m_opts.m_remap_ratio) return false;
    if (m_opts.m_remap_nlookup!=other.m_opts.m_remap_nlookup) return false;
    if (m_buckets!=other.m_buckets) return false;
    return true;
}

bool MappedTableBase::operator!=(const MappedTableBase &other) const {
    return !(*this==other);
}

uint_t MappedTableBase::nbucket_guess(uint_t nitem, double ratio) {
    return nitem / (2*ratio + 1);
}

uint_t MappedTableBase::nbucket() const {
    return m_buckets.size();
}

void MappedTableBase::clear_map() {
    for (auto &bucket: m_buckets) bucket.clear();
}

bool MappedTableBase::remap_due() const {
    return (m_nlookup_total >= m_opts.m_remap_nlookup) &&
        (double(m_nskip_total) / double(m_nlookup_total)) > m_opts.m_remap_ratio;
}

bool MappedTableBase::all_nonzero_records_mapped(const TableBase &source) const {
    // make a set out of all bucketed indices
    std::set<uint_t> set;
    for (const auto& bucket: m_buckets){
        for (const auto& item: bucket) set.insert(item);
    }
    // then go through table, if the row is non-zero, its index should be in the set.
    for (uint_t irow=0ul; irow<source.m_hwm; ++irow){
        bool in_set = set.find(irow)!=set.end();
        if (source.is_cleared(irow)==in_set) return false;
    }
    return true;
}