//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <forward_list>
#include <set>

#include <M7_lib/hdf5/Node.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/conf/Conf.h>

#include "Table.h"


struct MappedTableOptions {
    static constexpr uint_t c_default_nbucket_init = 100ul;
    double m_remap_ratio = conf::HashMapping::c_default_remap_ratio;
    uint_t m_remap_nlookup = conf::HashMapping::c_default_remap_nlookup;
    uint_t m_nbucket_init = c_default_nbucket_init;
    MappedTableOptions() = default;
    MappedTableOptions(const conf::HashMapping& opts) {
        m_remap_ratio = opts.c_default_remap_ratio;
        m_remap_nlookup = opts.m_remap_nlookup;
    }
};

struct MappedTableBase {
    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    mutable uint_t m_nskip_total = 0.0;
    mutable uint_t m_nlookup_total = 0.0;

    v_t<std::forward_list<uint_t>> m_buckets;

    MappedTableOptions m_mapping_opts;

    MappedTableBase(MappedTableOptions opts);

    MappedTableBase& operator=(const MappedTableBase& other);

    MappedTableBase(const MappedTableBase& other);

    bool operator==(const MappedTableBase& other) const;

    bool operator!=(const MappedTableBase& other) const;

    /**
     * @return
     *  size of the m_buckets vector
     */
    uint_t nbucket() const;
    /**
     * clear all buckets
     */
    void clear_map();
    /**
     * debugging only - checks that all nonzero records below the hwm of the source are mapped in the current m_buckets.
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
     * debugging only - checks that all nonzero records below the hwm of the source are mapped in the current m_buckets.
     * Defined here to reduce bloat of the templated class MappedTable
     * @param source
     *  base class cast of the mapped table
     * @return
     *  true if table passes verification
     */
    bool all_nonzero_records_mapped(const TableBase &source) const;

};

struct Lookup {
    const std::forward_list<uint_t>& m_bucket;
    std::forward_list<uint_t>::const_iterator m_prev;
    operator bool () const {
        return m_prev!=m_bucket.cend();
    }

    uint_t irow() const {
        if (!*this) return ~0ul;
        auto it = m_prev;
        ++it;
        return *it;
    }
};


template<typename row_t>
struct MappedTable : Table<row_t>, MappedTableBase {
    /**
     * won't compile unless the row defines a key_field_t;
     */
    typedef typename row_fields::Key<row_t>::type key_field_t;

    using Table<row_t>::to_string;
    using Table<row_t>::m_row;
    using TableBase::m_bw;

    row_t m_lookup_row;
    row_t m_insert_row;
    row_t m_erase_row;

    MappedTable(const row_t &row, MappedTableOptions opts) :
            Table<row_t>(row), MappedTableBase(opts),
            m_lookup_row(m_row), m_insert_row(m_row), m_erase_row(m_row){}

    MappedTable(const row_t &row) : MappedTable(row, MappedTableOptions()){}

    MappedTable& operator=(const MappedTable& other) {
        MappedTableBase::operator=(other);
        Table<row_t>::operator=(other);
        return *this;
    }

    MappedTable(const MappedTable<row_t> &other) : MappedTable(other.m_row, other.m_mapping_opts) {}

    bool operator==(const MappedTable<row_t> &other) const {
        if (!MappedTableBase::operator==(other)) return false;
        return TableBase::operator==(other);
    }

private:
    /**
     * attempt to find the element identified by key in the hash map.
     * @param key
     *  key to lookup
     * @return
     *  result object
     */
    Lookup lookup(const key_field_t &key, const v_t<std::forward_list<uint_t>>& buckets, const row_t& row) const {
        m_nlookup_total++;
        const auto nbucket = buckets.size();
        auto& bucket = buckets[key.hash() % nbucket];
        Lookup res = {bucket, bucket.cbefore_begin()};
        for (auto current = std::next(res.m_prev); current != bucket.cend(); current++) {
            row.jump(*current);
            if (row_fields::key(row) == key) return res;
            res.m_prev++;
            m_nskip_total++;
        }
        /*
         * key not found in map
         */
        static_cast<const Row&>(row).select_null();
        res.m_prev = bucket.cend();
        return res;
    }

public:
    Lookup lookup(const key_field_t& key, const row_t& row) const {
        return lookup(key, m_buckets, row);
    }

    const row_t& lookup(const key_field_t& key) const {
        lookup(key, m_lookup_row);
        return m_lookup_row;
    }

    row_t& lookup(const key_field_t& key) {
        lookup(key, m_lookup_row);
        return m_lookup_row;
    }

private:
    /**
     * same as the above function, but does not affect hash table statistics - useful for testing and ASSERTs
     * @param key
     *  key to lookup
     */
    Lookup uncounted_lookup(const key_field_t &key, const v_t<std::forward_list<uint_t>>& buckets, const row_t& row) {
        const auto nlookup_total = m_nlookup_total;
        const auto nskip_total = m_nskip_total;
        const auto res = lookup(key, buckets, row);
        m_nlookup_total = nlookup_total;
        m_nskip_total = nskip_total;
        return res;
    }

    Lookup uncounted_lookup(const key_field_t &key, const row_t& row) {
        return uncounted_lookup(key, m_buckets, row);
    }

public:
    void clear() override {
        TableBase::clear();
        clear_map();
    }

    void clear_row(row_t& row) {
        DEBUG_ASSERT_TRUE(Table<row_t>::associated(row), "row object not associated with this table");
        auto lookup = this->lookup(row_fields::key(row), row);
        erase(lookup);
    }

    /**
     * remove the record pointed to by row by:
     *  1. physically clearing by zeroing its record in the buffer,
     *  2. adding its index to the empty records stack
     *  3. removing its node in the associated bucket
     * @param lookup
     *  points to the bucket and element within it to be deleted
     */
    void erase(Lookup lookup) {
        DEBUG_ASSERT_TRUE(lookup, "cannot erase non-existent record");
        /*
         * do steps 1 and 2.
         */
        TableBase::free(lookup.irow());
        /*
         * for step 3 we must convert to non-const
         */
        auto& bucket = const_cast<std::forward_list<uint_t>&>(lookup.m_bucket);
        bucket.erase_after(lookup.m_prev);
    }

    void erase(const key_field_t &key) {
        erase(lookup(key, m_erase_row));
    }

    /**
     * insert into the hash table a new mapped element defined by the key arg
     * @param key
     *  hash table key object
     * @param row
     *  row object which points to the newly inserted record upon return of the method
     */
    void insert(const key_field_t &key, row_t& row) {
        DEBUG_ASSERT_TRUE(Table<row_t>::associated(row), "row object not associated with this table");
        DEBUG_ASSERT_FALSE(uncounted_lookup(key, row), "cannot insert when the key already exists in the table");
        auto irow = TableBase::get_free_row();
        auto &bucket = m_buckets[key.hash() % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        row.jump(irow);
        row.key_field() = key;
    }

    row_t& insert(const key_field_t &key) {
        insert(key, m_insert_row);
        return m_insert_row;
    }

    void insert(const row_t& src, row_t& row) {
        insert(row_fields::key(src), row);
        row.copy_in(src);
    }

    row_t& insert(const row_t& src) {
        insert(src, m_insert_row);
        return m_insert_row;
    }


    row_t& lookup_or_insert(const key_field_t& key) {
        lookup(key, m_lookup_row);
        if (m_lookup_row) return m_lookup_row;
        return insert(key);
    }

    row_t& lookup_or_insert(const row_t& src) {
        lookup(src, m_lookup_row);
        if (m_lookup_row) return m_lookup_row;
        return insert(src);
    }

private:
    /**
     * for an already-inserted row, map it to the hash table defined by the buckets arg
     * @param irow
     *  index of the row to be mapped
     * @param buckets
     *  vector of buckets into which the row indices are stored
     */
    void post_insert_buckets(uint_t irow, v_t<std::forward_list<uint_t>> &buckets) {
        DEBUG_ASSERT_FALSE(TableBase::m_is_freed_row[irow], "cannot map to a freed slot");
        auto& row = m_insert_row;
        row.jump(irow);
        const auto& key = row_fields::key(row);
        // row mustn't have already been added
        DEBUG_ASSERT_FALSE(uncounted_lookup(key, buckets, m_row),
                           "cannot map to a row which is already in the hash table");
        auto &bucket = buckets[key.hash() % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }
    /**
     * applies the general case by supplying the current m_buckets vector to the post_insert_buckets method
     * @param irow
     *  index of the row to be mapped
     */
    void post_insert(uint_t irow) {
        post_insert_buckets(irow, m_buckets);
    }

    /**
     * If a block of records was copied contiguously into a table with non-trivial post-insert obligations, a call to
     * post_insert is needed for each copied record
     * @param ibegin
     *  record index of beginning of the already-inserted data
     * @param iend
     *  record index of end of the already-inserted data
     */
    void post_insert_range(uint_t ibegin = 0ul, uint_t iend = ~0ul) {
        if (iend == ~0ul) iend = this->nrow_in_use();
        for (uint_t i = ibegin; i < iend; ++i) post_insert(i);
    }
public:
    /**
     * construct a new vector of buckets with a different size
     */
    void remap() {
        DEBUG_ASSERT_TRUE(all_nonzero_records_mapped(*this), "mapping is inconsistent with table row content");
        uint_t nbucket_new = nbucket() * skip_lookup_ratio() / m_mapping_opts.m_remap_ratio;
        // use the same expansion factor as for the Table buffer
        nbucket_new *= 1.0 + this->m_bw.get_expansion_factor();
        if (!TableBase::name().empty()) {
            logging::info_("remapping hash table for \"{}\"", TableBase::name());
            logging::info_("current ratio of skips to total lookups ({}) exceeds set limit ({})",
                           skip_lookup_ratio(), m_mapping_opts.m_remap_ratio);
            logging::info_("replacing current bucket vector of size {} with a new one of size {}",
                       nbucket(), nbucket_new);
        }
        v_t<std::forward_list<uint_t>> new_buckets(nbucket_new);
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

    void resize(uint_t nrec, double factor=-1.0) override {
        TableBase::resize(nrec, factor);
        attempt_remap();
    }

    void all_gatherv(const TableBase& src) override {
        clear();
        TableBase::all_gatherv(src);
        post_insert_range();
    }

    void gatherv(const TableBase& src, uint_t irank) override {
        clear();
        TableBase::gatherv(src);
        if (mpi::i_am(irank)) post_insert_range();
    }

    void load_fn(const hdf5::NodeReader& parent, const str_t& name, uint_t max_nitem_per_op,
                 strm_t field_names, bool part, bool this_rank) override {
        const auto key_name = m_row.key_field().m_name;
        auto it = std::find_if(field_names.cbegin(), field_names.cend(),
                               [&key_name](const strp_t& pair){return key_name == pair.first;});
        // the key field is not optional - it must be loaded
        if (it == field_names.cend()) field_names.insert({key_name, key_name});
        clear();
        Table<row_t>::load_fn(parent, name, max_nitem_per_op, field_names, part, this_rank);
        // run through all the loaded rows and add their key fields into the hash map
        for (uint_t irow = 0ul; irow < this->nrow_in_use(); ++irow) post_insert(irow);
    }
};

#endif //M7_MAPPEDTABLE_H
