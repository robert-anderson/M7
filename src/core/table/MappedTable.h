//
// Created by rja on 01/11/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include <forward_list>
#include <stack>
#include <src/core/parallel/RankAllocator.h>
#include "Table.h"
#include "src/core/field/TableField.h"
#include "src/core/table/BufferedTable.h"

struct LookupResult {
    std::forward_list<size_t> &m_bucket;
    std::forward_list<size_t>::iterator m_prev;

    LookupResult(std::forward_list<size_t> &bucket) :
            m_bucket(bucket), m_prev(bucket.before_begin()) {}

    operator bool() const {
        return m_prev != m_bucket.end();
    }

    size_t operator*() const {
        if (!*this) return ~0ul;//throw std::runtime_error("Cannot return stored index, hashmap entry not found!");
        return *std::next(m_prev);
    }
};

template<typename table_t, typename field_t, typename hash_fn=typename field_t::hash_fn>
struct MappedTable : table_t {
    static_assert(std::is_base_of<Table, table_t>::value, "Template arg must be derived from Table");

    using Table::push_back;
    using Table::get_free_row;
    using Table::m_row_size;
    using Table::m_row_dsize;
    using Table::m_bw;
    using Table::clear;
    using Table::dbegin;

    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 2.0;
    static constexpr size_t c_nbucket_min = 100;
    static_assert(std::is_base_of<NdFieldGroup<0ul>, field_t>::value, "Key field must be scalar");

    field_t &m_key_field;
    std::vector<std::forward_list<size_t>> m_buckets;

    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_total_skips = 0.0;
    size_t m_total_lookups = 0.0;
    double m_remap_ratio = c_default_remap_ratio;
    size_t m_remap_nlookup = c_default_remap_nlookup;

    template<typename ...Args>
    MappedTable(field_t &key_field, Args... args) :
            table_t(args...), m_key_field(key_field), m_buckets(c_nbucket_min) {}

    template<typename ...Args>
    MappedTable(size_t nbucket, field_t &key_field, Args... args) :
            table_t(args...), m_key_field(key_field), m_buckets(nbucket) {}

    defs::hash_t hash(typename field_t::view_t key) {
        return hash_fn()(key);
    }

    typename field_t::view_t get_key(const size_t &irow) {
        return m_key_field(irow);
    }

    LookupResult operator[](const typename field_t::view_t key) {
        m_total_lookups++;
        LookupResult res(m_buckets[hash(key) % m_buckets.size()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            if (get_key(*next) == key) return res;
            res.m_prev++;
            m_total_skips++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    using typename Table::cb_list_t;

    void erase_rows(const defs::inds &irows, const cb_list_t &callbacks) override {
        for (auto irow : irows) {
            auto lookup = (*this)[m_key_field(irow)];
            erase(lookup);
        }
    }

    void insert_rows(const Buffer &recv, const cb_list_t &callbacks) override {
        const auto nrow = recv.dsize() / m_row_dsize;
        for (size_t irow = 0; irow < nrow; ++irow) {
            auto iinsert = get_free_row();
            std::memcpy(dbegin(iinsert), recv.dbegin() + irow * m_row_dsize, m_row_size);
            post_insert(iinsert);
        }
    }

    void erase(LookupResult result) {
        clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }


    size_t insert(const typename field_t::view_t key) {
        ASSERT(!(this->operator[](key)))
        auto irow = get_free_row();
        auto &bucket = m_buckets[hash(key) % m_buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_key_field(irow) = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!Table::is_cleared(irow))
        auto key = m_key_field(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[key]);
        auto &bucket = buckets[hash(key) % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) {
        post_insert(irow, m_buckets);
    }

    void print_map() const {
        size_t nkey = 0ul;
        for (size_t ibucket = 0ul; ibucket < m_buckets.size(); ++ibucket) {
            auto bucket = m_buckets[ibucket];
            if (bucket.empty()) continue;
            std::cout << "Bucket " << std::to_string(ibucket) << std::endl;
            for (auto irow: bucket) {
                std::cout << "\t" << m_key_field(irow).to_string() << " => " << irow << std::endl;
                nkey++;
            }
        }
        std::cout << "Number of keys: " << nkey << std::endl;
    }

    void remap_if_required() {
        if (m_total_lookups < m_remap_nlookup) return;
        if (double(m_total_skips) / double(m_total_lookups) > m_remap_ratio) return;
        // use the same expansion factor as for the Table buffer
        const size_t nbucket_new = m_buckets.size() * m_bw.expansion_factor();
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

    struct DynamicRowSet : RankAllocator<field_t>::Dynamic {
        static_assert(std::is_base_of<Table, table_t>::value, "Template arg must be derived from Table");
        typedef RankAllocator<field_t> ra_t;
        typedef MappedTable<table_t, field_t, hash_fn> mt_t;
        /*
         * the mapped table which stores the definitive row values
         */
        const mt_t &m_source;
        /*
         * the unmapped table which loads copies of rows between m_source (arbitrary order, non-
         * contiguous) and m_all (contiguous). lc = "local, contiguous"
         */
        BufferedTable<table_t> m_lc;
        /*
         * the unmapped table which holds copies of all mapped rows. these copies are refreshed
         * with a call to refresh method. This table is intended for reading rows from all MPI ranks,
         * as such there is no machanism for committing changes to m_source from m_all. It is
         * to be treated as a read-only copy. ac = "all, contiguous"
         */
        BufferedTable<table_t> m_ac;
        /*
         * map row index in source -> row index in m_local.
         *
         * We could have achieved this mapping by extending the table_t type with an additional
         * field of type fields::Number<size_t>, then extending MappedTable with the resultant
         * type as a template argument. This has the downside that rows cannot be copied wholesale
         * from m_source into/out of such a table, since the rows would be longer in the extended
         * table_t. Better to use the stl-provided map to point into a table of the same format.
         */
        std::map<size_t, size_t> m_map;
        /*
         * once gathered, the major ordering of m_all will be MPI rank index ascending. m_counts stores
         * the number of rows tracked on each rank and m_displs stores the row index of m_all where
         * each rank's rows begin.
         */
        defs::inds m_counts;
        defs::inds m_displs;
        /*
         * keep track of which ranks have any rows
         */
        bool m_have_any_rows;
        defs::inds m_ranks_with_any_rows;

    public:
        DynamicRowSet(ra_t &ra, const mt_t &mt) :
                ra_t::Dynamic(ra),
                m_source(mt),
                m_lc("Rank-local dynamic row set rows", static_cast<const table_t &>(mt)),
                m_ac("All dynamic row set rows", static_cast<const table_t &>(mt)),
                m_counts(mpi::nrank(), 0ul),
                m_displs(mpi::nrank(), 0ul) {
            m_ranks_with_any_rows.reserve(mpi::nrank());
        }

        size_t nrow() const {
            return m_map.size();
        }

        void add(size_t irow) {
            m_map[irow] = nrow();
        }

    private:
        void populate_local() {
            static_cast<Table &>(m_lc).clear();
            m_lc.push_back(nrow());
            /*
             * local row indices must be updated in case we sent a block
             */
            size_t irow_local = 0ul;
            for (auto &pair : m_map) {
                pair.second = irow_local++;
                static_cast<Table &>(m_lc).copy_row_in(m_source, pair.first, pair.second);
            }
        }

        void gatherv() {
            m_ac.clear();
            mpi::all_gather(nrow(), m_counts);
            mpi::counts_to_displs_consec(m_counts, m_displs);
            auto nrow = m_displs.back() + m_counts.back();
            m_ac.push_back(nrow);
            /*
             * convert from units of rows to datawords...
             */
            for (auto &i : m_counts) i *= m_source.m_row_dsize;
            for (auto &i : m_displs) i *= m_source.m_row_dsize;
            mpi::all_gatherv(m_lc.dbegin(), m_counts[mpi::irank()], m_ac.dbegin(), m_counts, m_displs);
            /*
             * ... and back again
             */
            for (auto &i : m_counts) i /= m_source.m_row_dsize;
            for (auto &i : m_displs) i /= m_source.m_row_dsize;
            m_have_any_rows = m_counts[mpi::irank()];
            m_ranks_with_any_rows.clear();
            for (auto &i : m_counts) {
                if (i) m_ranks_with_any_rows.push_back(mpi::irank());
            }
        }

    public:
        virtual void post_update() {}
        void update() {
            populate_local();
            gatherv();
            post_update();
        }

        void on_row_send_(size_t irow) override {
            m_map.erase(irow);
        }

        void on_row_recv_(size_t irow) override {
            ASSERT(m_map.find(irow) == m_map.end());
            m_map[irow] = nrow();
        }

    };

    DynamicRowSet dynamic_row_set(RankAllocator<field_t> &ra) {
        return DynamicRowSet(ra, *this);
    }


    struct DynamicRow : public DynamicRowSet {
        using typename DynamicRowSet::ra_t;
        using typename DynamicRowSet::mt_t;
        using ra_t::Dynamic::m_ra;
        using DynamicRowSet::m_source;
        using DynamicRowSet::m_map;
        using DynamicRowSet::add;
        using DynamicRowSet::nrow;
        using DynamicRowSet::update;
        using DynamicRowSet::m_have_any_rows;
        using DynamicRowSet::m_ranks_with_any_rows;

        DynamicRow(ra_t &ra, const mt_t &mt, Table::Loc loc) :
                DynamicRowSet(ra, mt) {
            if (loc.is_mine()) add(loc.m_irow);
            update();
            if (m_ranks_with_any_rows.size()!=1) mpi::stop_all("Only one rank should have a row");
        }

        void change(Table::Loc loc) {
            if (loc.is_mine()) {
                ASSERT(m_ra.get_rank(m_source.m_key_field(loc.m_irow))==mpi::irank());
                m_map.clear();
                m_map = {loc.m_irow, 0ul};
            }
            update();
        }

        bool is_mine() const {
            ASSERT(m_have_any_rows);
            ASSERT(m_ranks_with_any_rows.size()==1);
            return DynamicRowSet::nrow();
        }
    };

};


#endif //M7_MAPPEDTABLE_H
