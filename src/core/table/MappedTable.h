//
// Created by rja on 10/02/2021.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <forward_list>
#include <set>
#include <io/HDF5Wrapper.h>
#include "Table.h"
#include "field/Fields.h"

struct LookupResult {
    std::forward_list<size_t> &m_bucket;
    std::forward_list<size_t>::iterator m_prev;

    LookupResult(std::forward_list<size_t> &bucket) :
            m_bucket(bucket), m_prev(bucket.before_begin()) {}

    operator bool() const {
        return m_prev != m_bucket.end();
    }

    size_t operator*() const {
        if (!*this) return ~0ul;
        return *std::next(m_prev);
    }
};


struct MappedTableBase {
    static constexpr size_t c_default_nbucket = 100;
    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 500;
    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_nskip_total = 0.0;
    size_t m_nlookup_total = 0.0;

    std::vector<std::forward_list<size_t>> m_buckets;
    const size_t m_remap_nlookup;
    const double m_remap_ratio;

    MappedTableBase(size_t nbucket = 0ul, size_t remap_nlookup = 0ul, double remap_ratio = 0.0);

    /**
     * Assuming uniform frequency of access, the average number of skips per lookup in a bucket containing n items is
     *  (0 + 1 + ... + n-1) / n = (n-1)/2
     * So assuming uniform distribution of items in buckets, the expected number of skips per lookup is
     *  (nitem/nbucket - 1) / 2
     * this is the remap ratio, so invert this estimate to obtain the required number of buckets
     * @param nitem
     *  number of items (nonzero table rows) on which to base the capacity estimation
     * @param ratio
     *  the maximum acceptable ratio of skips to lookups
     * @return
     *  estimated number of buckets required to support the mapping of nitem items with the given skip ratio
     */
    static size_t nbucket_guess(size_t nitem, double ratio);
    /**
     * @return
     *  size of the m_buckets vector
     */
    size_t nbucket() const;
    /**
     * clear all buckets
     */
    void clear_map();
    /**
     * debugging only - checks that all nonzero rows below the hwm of the source are mapped in the current m_buckets.
     * Defined here to reduce bloat of the templated class MappedTable
     * @param source
     *  base class cast of the mapped table
     * @return
     *  true if table passes verification
     */
    double skip_lookup_ratio() const {
        return double(m_nskip_total)/double(m_nlookup_total);
    }
    /**
     * @return
     *  true if enough lookups have been attempted and the ratio of skips/lookups is worse than m_remap_ratio
     */
    bool remap_due() const;
    /**
     * debugging only - checks that all nonzero rows below the hwm of the source are mapped in the current m_buckets.
     * Defined here to reduce bloat of the templated class MappedTable
     * @param source
     *  base class cast of the mapped table
     * @return
     *  true if table passes verification
     */
    bool all_nonzero_rows_mapped(const TableBase &source) const;

};


template<typename row_t>
struct MappedTable : Table<row_t>, MappedTableBase {
    /*
     * won't compile unless the row defines a key_field_t;
     */
    typedef typename KeyField<row_t>::type key_field_t;
    row_t m_lookup_row;
    row_t m_insert_row;

    using Table<row_t>::m_row;
    using Table<row_t>::to_string;
    using TableBase::m_hwm;
    using TableBase::m_bw;

    MappedTable(const row_t &row, size_t nbucket = 0ul, size_t remap_nlookup = 0ul, double remap_ratio = 0.0) :
            Table<row_t>(row), MappedTableBase(nbucket, remap_nlookup, remap_ratio),
            m_lookup_row(m_row), m_insert_row(m_row) {}

    MappedTable(const MappedTable<row_t> &other) : MappedTable(other.m_row, other.m_buckets.size(),
                                                               other.m_remap_nlookup, other.m_remap_ratio) {}
    /**
     * attempt to find the element identified by key in the hash map.
     * @param key
     *  key to lookup
     */
    LookupResult lookup(const key_field_t &key, std::vector<std::forward_list<size_t>>& buckets) {
        m_nlookup_total++;
        auto nbucket = buckets.size();
        LookupResult res(buckets[key.hash() % nbucket]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            m_lookup_row.jump(*next);
            if (KeyField<row_t>::get(m_lookup_row) == key) return res;
            res.m_prev++;
            m_nskip_total++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    LookupResult operator[](const key_field_t &key) {
        return lookup(key, m_buckets);
    }

private:
    /**
     * same as the above function, but does not affect hash table statistics - useful for testing and ASSERTs
     * @param key
     *  key to lookup
     */
    LookupResult uncounted_lookup(const key_field_t &key) {
        auto nlookup_total = m_nlookup_total;
        auto nskip_total = m_nskip_total;
        auto lookup = (*this)[key];
        m_nlookup_total = nlookup_total;
        m_nskip_total = nskip_total;
        return lookup;
    }

public:
    void clear() override {
        TableBase::clear();
        clear_map();
    }

    void clear(const size_t &irow) override {
        m_row.jump(irow);
        auto lookup = (*this)[KeyField<row_t>::get(m_row)];
        erase(lookup);
    }

    /**
     * remove the row pointed to by lookup by:
     *  1. physically clearing the row by zeroing its slice of the buffer,
     *  2. adding its index to the free rows stack
     *  3. removing its node in the associated bucket
     * @param lookup
     *  points to the bucket element to be deleted
     */
    void erase(LookupResult &lookup) {
        TableBase::clear(*lookup);
        lookup.m_bucket.erase_after(lookup.m_prev);
        // put into "not found" state:
        lookup.m_prev = lookup.m_bucket.end();
    }
    /**
     * insert into the hash table a new mapped element defined by the key arg
     * @param key
     *  hash table key object
     * @return
     *  row index of the key
     */
    size_t insert(const key_field_t &key) {
        DEBUG_ASSERT_FALSE(uncounted_lookup(key), "cannot insert when the key already exists in the table");
        auto irow = TableBase::get_free_row();
        auto &bucket = m_buckets[key.hash() % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_row.jump(irow);
        m_row.key_field() = key;
        return irow;
    }
    /**
     * for an already-inserted row, map it to the hash table defined by the buckets arg
     * @param irow
     *  index of the row to be mapped
     * @param buckets
     *  vector of buckets into which the row indices are stored
     */
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        DEBUG_ASSERT_FALSE(TableBase::is_cleared(irow), "cannot map to a cleared row");
        m_insert_row.jump(irow);
        // row mustn't have already been added
        DEBUG_ASSERT_FALSE(lookup(KeyField<row_t>::get(m_insert_row), buckets),
                           "cannot map to a row which is already in the hash table");
        auto &bucket = buckets[KeyField<row_t>::get(m_insert_row).hash() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }
    /**
     * applies the general case by supplying the current m_buckets vector to the post_insert_buckets method
     * @param irow
     *  index of the row to be mapped
     */
    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }
    /**
     * construct a new vector of buckets with a different size
     */
    void remap() {
        DEBUG_ASSERT_TRUE(all_nonzero_rows_mapped(*this), "mapping is inconsistent with table row content");
        size_t nbucket_new = nbucket() * skip_lookup_ratio()/m_remap_ratio;
        // use the same expansion factor as for the Table buffer
        nbucket_new *= 1.0 + this->m_bw.get_expansion_factor();
        if (!TableBase::name().empty()) {
            log::info_("remapping hash table for \"{}\"", TableBase::name());
            log::info_("current ratio of skips to total lookups ({}) exceeds set limit ({})",
                       skip_lookup_ratio(), m_remap_ratio);
            log::info_("replacing current bucket vector of size {} with a new one of size {}",
                       nbucket(), nbucket_new);
        }
        std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
        for (const auto &old_bucket : m_buckets) {
            for (const auto irow : old_bucket) {
                post_insert_buckets(irow, new_buckets);
            }
        }
        // destroy current m_buckets and put new_buckets in its place
        m_buckets = std::move(new_buckets);
        DEBUG_ASSERT_EQ(nbucket(), nbucket_new, "error in bucket vector update");
    }
    /**
     * if a remap is due (enough lookups + bad enough skips/lookups ratio), the remap will be executed and those
     * counters reset
     */
    void attempt_remap() {
        if (remap_due()) {
            remap();
            // reset counters
            m_nskip_total = 0ul;
            m_nlookup_total = 0ul;
        }
    }

    void load(hdf5::GroupReader &parent, std::string name) override {
        RowHdf5Reader<row_t> row_reader(m_row, parent, name);
        size_t iitem = 0ul;
        clear();
        TableBase::push_back(row_reader.m_nitem);
        for (row_reader.restart(); row_reader.in_range(); row_reader.step()) {
            row_reader.read(iitem++);
            post_insert(row_reader.index());
        }
    }

};

#endif //M7_MAPPEDTABLE_H
