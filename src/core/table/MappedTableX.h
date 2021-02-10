//
// Created by rja on 07/02/2021.
//

#ifndef M7_MAPPEDTABLEX_H
#define M7_MAPPEDTABLEX_H

#include <forward_list>
#include "src/core/field/FieldX.h"
#include "TableBaseX.h"

struct LookupResultX {
    std::forward_list<size_t> &m_bucket;
    std::forward_list<size_t>::iterator m_prev;

    LookupResultX(std::forward_list<size_t> &bucket) :
            m_bucket(bucket), m_prev(bucket.before_begin()) {}

    operator bool() const {
        return m_prev != m_bucket.end();
    }

    size_t operator*() const {
        if (!*this) return ~0ul;
        return *std::next(m_prev);
    }
};


template<typename row_t>
struct MappedTableX : TableX<row_t> {
    static_assert(std::is_base_of<RowX, row_t>::value, "Template arg must be derived from Row");
    using TableX<row_t>::m_row;
    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 2.0;
    static constexpr size_t c_nbucket_min = 100;
    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_total_skips = 0.0;
    size_t m_total_lookups = 0.0;
    double m_remap_ratio = c_default_remap_ratio;
    size_t m_remap_nlookup = c_default_remap_nlookup;

    std::vector<std::forward_list<size_t>> m_buckets;

    MappedTableX(row_t row, size_t nbucket):
            TableX<row_t>(row), m_buckets(nbucket){
        MPI_ASSERT(m_row.m_key_field_inds.size(), "A Table can't be mapped without key fields");
    }

    /**
     * @param nrow_max
     * expected maximum number of rows
     * @param mean_skips_per_bucket
     * the maximum acceptable value for the average number
     * @return
     */
    static size_t nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket) {
        // number of skips is on average half the average number of rows per bin.
        return std::max(c_nbucket_min, size_t(nrow_max / double(2 * mean_skips_per_bucket)));
    }

    LookupResultX operator[](const RowX& key) {
        m_total_lookups++;
        LookupResultX res(m_buckets[key.hash_keys() % m_buckets.size()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            m_row.select(*next);
            if (m_row.equal_keys(key)) return res;
            res.m_prev++;
            m_total_skips++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    void erase_rows(const defs::inds &irows) override {
        for (auto irow : irows) {
            m_row.select(irow);
            auto lookup = (*this)[m_row];
            erase(lookup);
        }
    }

    void erase(LookupResultX& result) {
        TableBaseX::clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }


    size_t insert(const RowX& row) {
        ASSERT(!(this->operator[](row)))
        auto irow = TableBaseX::get_free_row();
        auto &bucket = m_buckets[row.hash_keys() % m_buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_row.select(irow);
        m_row.copy_keys(row);
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!TableBaseX::is_cleared(irow))
        m_row.select(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[m_row]);
        auto &bucket = buckets[m_row.hash_keys() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }

    void remap_if_required() {
        if (m_total_lookups < m_remap_nlookup) return;
        if (double(m_total_skips) / double(m_total_lookups) > m_remap_ratio) return;
        // use the same expansion factor as for the Table buffer
        const size_t nbucket_new = m_buckets.size() * TableBaseX::m_bw.expansion_factor();
        std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
        for (const auto &old_bucket : m_buckets) {
            for (const auto irow : old_bucket) {
                post_insert(irow, new_buckets);
            }
        }
        m_buckets = std::move(new_buckets);
        // reset counters
        m_total_skips = 0ul;
        m_total_lookups = 0ul;
    }

};


#endif //M7_MAPPEDTABLEX_H
