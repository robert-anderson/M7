//
// Created by rja on 10/02/2021.
//

#ifndef M7_MAPPEDTABLEZ_H
#define M7_MAPPEDTABLEZ_H

#include <forward_list>
#include "TableZ.h"

struct LookupResultZ {
    std::forward_list<size_t> &m_bucket;
    std::forward_list<size_t>::iterator m_prev;

    LookupResultZ(std::forward_list<size_t> &bucket) :
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
struct MappedTableZ : TableZ<row_t> {
    /*
     * won't compile unless the row defines a key_field_t;
     */
    typedef decltype(row_t::m_key_field) key_field_t;

    using TableZ<row_t>::m_row;

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


    MappedTableZ(const row_t& row, size_t nbucket) : TableZ<row_t>(row), m_buckets(nbucket){}

    size_t nbucket() const {
        return nbucket();
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

    LookupResultZ operator[](const key_field_t& key) {
        m_total_lookups++;
        LookupResultZ res(m_buckets[key.hash() % nbucket()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            m_row.jump(*next);
            if (m_row.m_key_field.equals(key)) return res;
            res.m_prev++;
            m_total_skips++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    void erase_rows(const defs::inds &irows) override {
        for (auto irow : irows) {
            m_row.jump(irow);
            auto lookup = (*this)[m_row.m_key_field];
            erase(lookup);
        }
    }

    void erase(LookupResultZ& result) {
        TableBaseZ::clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }


    size_t insert(const key_field_t& key) {
        ASSERT(!(this->operator[](key)));
        auto irow = TableBaseZ::get_free_row();
        auto &bucket = m_buckets[key.hash() % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_row.select(irow);
        m_row.m_key_field = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!TableBaseZ::is_cleared(irow))
        m_row.jump(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[m_row.m_key_field]);
        auto &bucket = buckets[m_row.m_key_field.hash() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }

    void remap_if_required() {
        if (m_total_lookups < m_remap_nlookup) return;
        if (double(m_total_skips) / double(m_total_lookups) > m_remap_ratio) return;
        // use the same expansion factor as for the Table buffer
        const size_t nbucket_new = nbucket() * TableBaseZ::m_bw.expansion_factor();
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


#endif //M7_MAPPEDTABLEZ_H
