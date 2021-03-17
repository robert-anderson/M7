//
// Created by rja on 10/02/2021.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <forward_list>
#include <src/core/io/HDF5Wrapper.h>
#include "Table.h"
#include "src/core/field/Fields.h"

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
    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 500;
    static constexpr size_t c_nbucket_min = 100;
    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_ntotal_skip = 0.0;
    size_t m_ntotal_lookup = 0.0;
    double m_remap_ratio = c_default_remap_ratio;
    size_t m_remap_nlookup = c_default_remap_nlookup;

    std::vector<std::forward_list<size_t>> m_buckets;

    MappedTableBase(size_t nbucket);

    /**
     * @param nrow_max
     * expected maximum number of rows
     * @param mean_skips_per_bucket
     * the maximum acceptable value for the average number
     * @return
     */
    static size_t nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket);

    size_t nbucket() const;

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

    MappedTable(const row_t &row, size_t nbucket) : Table<row_t>(row), MappedTableBase(nbucket),
        m_lookup_row(m_row), m_insert_row(m_row) {}

    MappedTable(const MappedTable<row_t>& other): MappedTable(other.m_row, other.m_buckets.size()){}

    LookupResult operator[](const key_field_t &key) {
        m_ntotal_lookup++;
        LookupResult res(m_buckets[key.hash() % nbucket()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            m_lookup_row.jump(*next);
            if (KeyField<row_t>::get(m_lookup_row) == key) return res;
            res.m_prev++;
            m_ntotal_skip++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    void erase_rows(const defs::inds &irows) override {
        for (auto irow : irows) {
            m_row.jump(irow);
            auto lookup = (*this)[KeyField<row_t>::get(m_row)];
            erase(lookup);
        }
    }

    void erase(LookupResult &result) {
        TableBase::clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }


    size_t insert(const key_field_t &key) {
        ASSERT(!(this->operator[](key)));
        auto irow = TableBase::get_free_row();
        auto &bucket = m_buckets[key.hash() % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_row.jump(irow);
        m_row.key_field() = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!TableBase::is_cleared(irow))
        m_insert_row.jump(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[KeyField<row_t>::get(m_insert_row)]);
        auto &bucket = buckets[KeyField<row_t>::get(m_insert_row).hash() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }

    void remap_if_required() {
        return;
        if (m_ntotal_lookup > m_remap_nlookup && (double(m_ntotal_skip) / double(m_ntotal_lookup)) > m_remap_ratio) {
            // use the same expansion factor as for the Table buffer
            const size_t nbucket_new = nbucket() * TableBase::m_bw.expansion_factor();
            std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
            for (const auto &old_bucket : m_buckets) {
                for (const auto irow : old_bucket) {
                    post_insert_buckets(irow, new_buckets);
                }
            }
            m_buckets = std::move(new_buckets);
        }
        // reset counters
        m_ntotal_skip = 0ul;
        m_ntotal_lookup = 0ul;
    }

    void save(hdf5::Group& parent, std::string name) {
        //hdf5::NdListWriter<>
    }

};

#if 0
template<typename row_t, typename mapped_row_t>
struct IndirectMappedTable : Table<row_t>, MappedTableBase {
    /*
     * won't compile unless the row_t and mapped_row_t defines a key_field_t;
     * and the row_t is an integral field
     *
     * to get access to hashable key (mapped_row_t::key_field()), the integer
     * number field row_t::key_field() must be used to extract the index from the
     * parent table
     *
     */
    template<typename T>
    void index_type_check(const fields::Number<T>& field){
        static_assert(std::is_integral<T>::value, "indexing type must be integral");
    }

    typedef typename KeyField<mapped_row_t>::type key_field_t;
    typedef typename KeyField<row_t>::type index_field_t;

    const Table<mapped_row_t>& m_mapped_table;
    mapped_row_t m_mapped_lookup_row;
    row_t m_lookup_row;
    mapped_row_t m_mapped_insert_row;
    row_t m_insert_row;

    using Table<row_t>::m_row;
    using Table<row_t>::to_string;
    using TableBase::m_hwm;
    using TableBase::m_bw;

    IndirectMappedTable(const Table<mapped_row_t>& mapped_table, const row_t &row, size_t nbucket) :
            Table<row_t>(row), MappedTableBase(nbucket),
            m_mapped_table(mapped_table), m_mapped_lookup_row(m_mapped_table.m_row),
            m_mapped_insert_row(m_mapped_table.m_row),
            m_lookup_row(m_row), m_insert_row(m_row) {
        index_type_check(KeyField<row_t>::get(m_row));
    }

    IndirectMappedTable(const IndirectMappedTable<row_t, mapped_row_t>& other):
            IndirectMappedTable(other.m_row, other.m_buckets.size()){}

    LookupResult operator[](const key_field_t &key) {
        m_ntotal_lookup++;
        LookupResult res(m_buckets[key.hash() % nbucket()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            m_lookup_row.jump(*next);
            m_mapped_lookup_row.jump(KeyField<row_t>::get(m_lookup_row));
            if (KeyField<row_t>::get(m_mapped_lookup_row) == key) return res;
            res.m_prev++;
            m_ntotal_skip++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    /*
    void erase_rows(const defs::inds &irows) override {
        for (auto irow : irows) {
            m_row.jump(irow);
            auto lookup = (*this)[KeyField<row_t>::get(m_row)];
            erase(lookup);
        }
    }

    void erase(LookupResult &result) {
        TableBase::clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }
     */


    size_t insert(const key_field_t &key) {
        ASSERT(!(this->operator[](key)));
        auto irow = TableBase::get_free_row();
        auto &bucket = m_buckets[key.hash() % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_insert_row.jump(irow);
        m_mapped_insert_row.jump(KeyField<row_t>::get(m_insert_row));
        KeyField<row_t>::get(m_mapped_insert_row) = key;
        return irow;
    }

    /*
    // for situations in which the row has already been copied-to
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!TableBase::is_cleared(irow))
        m_insert_row.jump(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[KeyField<row_t>::get(m_insert_row)]);
        auto &bucket = buckets[KeyField<row_t>::get(m_insert_row).hash() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }

    void remap_if_required() {
        return;
        if (m_ntotal_lookup > m_remap_nlookup && (double(m_ntotal_skip) / double(m_ntotal_lookup)) > m_remap_ratio) {
            // use the same expansion factor as for the Table buffer
            const size_t nbucket_new = nbucket() * TableBase::m_bw.expansion_factor();
            std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
            for (const auto &old_bucket : m_buckets) {
                for (const auto irow : old_bucket) {
                    post_insert_buckets(irow, new_buckets);
                }
            }
            m_buckets = std::move(new_buckets);
        }
        // reset counters
        m_ntotal_skip = 0ul;
        m_ntotal_lookup = 0ul;
    }
     */

};

#endif //M7_MAPPEDTABLE_H
#endif //M7_MAPPEDTABLE_H
