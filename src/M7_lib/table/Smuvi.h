//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_SMUVI_H
#define M7_SMUVI_H

#include "MappedTable.h"
#include "BufferedTable.h"
#include "M7_lib/communication/SendRecv.h"
#include "M7_lib/communication/Distribution.h"


/**
 * Shared Memory Unordered map to Vectors of Items "SMUVI"
 * The SMUVI is a type of parallel hash map from a key domain to a variable-length vector of indices which has fewer
 * requirements than a parallelised generalization of a map from type T to a arbitrarily sized array of type U elements,
 * expressed using STL containers as std::unordered_map<T, std::vector<unsigned long>>.
 *
 * Primary among which is that the map is never simultaneously modified and read; it is at any given time in one of two
 * states, "insert" or "access".
 *
 * The implemented data structure need not support some operations commonly required by serial and concurrent hash maps:
 *  - Reassignment of existing elements. Newly inserted scalar values are appended to the existing list of values
 *    corresponding to the key
 *  - Deletion. No key or element of a value array may be removed after it has been added
 *  - Dynamic resizing of the "access" key array.
 *
 * The SMUVI behaves as a hash table which can be inserted to in parallel, then __collated__ with minimal
 * synchronisation overhead, to then be accessed for reading in constant time among all ranks in a communicator with a
 * common shared memory window.
 *
 * The working of this data structure is perhaps best explained by an example:
 *
 * Suppose we have 4 MPI ranks in 2 shared memory realms initially in the insert state
 * rank 0 (shmem 0, irank_shmem 0)
 * rank 1 (shmem 0, irank_shmem 1)
 * rank 2 (shmem 1, irank_shmem 0)
 * rank 3 (shmem 1, irank_shmem 1)
 *
 * Suppose we have the following insertion workloads on each rank:
 *  rank 0: [(C, 1), (B, 2), (A, 4), (B, 3), (A, 2) ]
 *  rank 1: [(A, 3), (C, 0), (A, 0), (D, 4) ]
 *  rank 2: [(A, 5), (D, 0), (B, 1), (C, 4) ]
 *  rank 3: [(B, 5), (D, 1), (D, 3), (C, 2) ]
 *
 * The final data accessible data structure must be equivalent to a copy of the following map:
 *  {
 *      A : [0, 2, 3, 4, 5 ],
 *      B : [1, 2, 3, 5 ],
 *      C : [0, 1, 2, 4 ],
 *      D : [0, 1, 3, 4 ]
 *  }
 *
 * Each key must be allotted a destination rank number in each shared memory region, suppose we have
 *  A -> [1, 3 ], B -> [0, 2 ], C -> [1, 2], D -> [0, 2]
 * these are determined by Distribution::one_irank_in_each_shmem_region
 *
 * Now consider the insertion procedure on a single process, say rank 0.
 * The inserter object is a communicating pair with nrank(=4) send tables and one receiving table
 * Initially the inserter send tables are empty:
 *  inserter send = [[], [], [], []]
 * then after the first key-index pair in the rank 0 workload - (C, 1) - they become
 *  inserter send = [[], [(C, 1)], [(C, 1)], []]
 * since those are the rank indices associated with key C.
 * repeating this process for the entire insertion workload of rank 0 we have
 *  inserter send = [[(B, 2), (B, 3)], [(C, 1), (A, 4), (A, 2)], [(C, 1), (B, 2), (B, 3)], [(A, 4), (A, 2)]]
 *
 * and on the other ranks:
 * rank 1:
 *  inserter send = [[(D, 4)], [(A, 3), (C, 0), (A, 0)], [(C, 0), (D, 4)], [(A, 3), (A, 0)]]
 * rank 2:
 *  inserter send = [[(D, 0), (B, 1)], [(A, 5), (C, 4)], [(D, 0), (B, 1), (C, 4)], [(A, 5)]]
 * rank 3:
 *  inserter send = [[(B, 5), (D, 1), (D, 3)], [(C, 2)], [(B, 5), (D, 1), (D, 3), (C, 2)], []]
 *
 * once the insertion is done, the collate method is called. The first process of which is to communicate the inserted
 * data to the receiving tables. This yields the following receiving tables
 * rank 0:
 *  inserter recv = [(B, 2), (B, 3), (D, 4), (D, 0), (B, 1), (B, 5), (D, 1), (D, 3)]
 * rank 1:
 *  inserter recv = [(C, 1), (A, 4), (A, 2), (A, 3), (C, 0), (A, 0), (A, 5), (C, 4), (C, 2)]
 * rank 2:
 *  inserter recv = [(C, 1), (B, 2), (B, 3), (C, 0), (D, 4), (D, 0), (B, 1), (C, 4), (B, 5), (D, 1), (D, 3), (C, 2)]
 * rank 3:
 *  inserter recv = [(A, 4), (A, 2), (A, 3), (A, 0), (A, 5)]
 *
 * Notice that each shared memory region (0, 1) and (2, 3) contains all key-index pairs between its constituent ranks.
 * Now each rank creates a key-(sorted indices) hash table *in shared memory* so that the rank can write without data
 * races, but once the collate step is complete, all ranks in the shared memory region can read from the created hash tables
 *
 * shmem 0:
 *  accessors = [
 *      { B: [1, 2, 3, 5 ], D : [0, 1, 3, 4 ] },
 *      { A: [0, 2, 3, 4, 5 ], C: [0, 1, 2, 4 ] }
 *  ]
 * shmem 1:
 *  accessors = [
 *      { C: [0, 1, 2, 4 ], B: [1, 2, 3, 5 ], D : [0, 1, 3, 4 ] },
 *      { A: [0, 2, 3, 4, 5 ] }
 *  ]
 *
 * at this point, the collate step is complete.
 *
 * After collation, the entries corresponding to each accessor table are stored in a flat shared memory array. Accesses
 * return a pointer to the beginning of the entries to which the key corresponds
 */


template<typename key_t, typename value_t>
class Smuvi {

    /**
     * Row type for the inserter tables
     */
    struct InsertRow : Row {
        key_t m_key;
        value_t m_value;

        InsertRow(const key_t& key, const value_t& value): m_key(key), m_value(value){
            static_cast<FieldBase&>(m_key).add_to_row(this);
            static_cast<FieldBase&>(m_value).add_to_row(this);
        }

        key_t &key_field() {
            return m_key;
        };
    };
    /**
     * Row type for the values tables
     */
    struct ValueRow : Row {
        value_t m_value;

        explicit ValueRow(const value_t& value): m_value(value){
            static_cast<FieldBase&>(m_value).add_to_row(this);
        }
    };
    /**
     * Row type for the accessor tables
     */
    struct AccessRow : Row {
        key_t m_key;
        field::Number<uint_t> m_value_count;
        field::Number<uint_t> m_value_displ;

        explicit AccessRow(const key_t& key): m_key(key), m_value_count(this), m_value_displ(this){
            static_cast<FieldBase&>(m_key).add_to_row(this);
        }

        key_t &key_field() {
            return m_key;
        };
    };
    /**
     * send and receive tables for the key-value pair insertion
     */
    send_recv::BasicSend<InsertRow> m_inserter;
    uintv_t m_irank_world_to_iaccess_table;
    /**
     * one shared-memory MappedTable for each rank in the shared memory region
     */
    v_t<buffered::MappedTable<AccessRow>> m_access_tables;
    v_t<AccessRow> m_access_foreach_rows;
    /**
     * a single flat storage for all the values associated with all keys in the shared memory region
     */
    v_t<buffered::Table<ValueRow>> m_values_tables;
    v_t<ValueRow> m_values_lookup_rows;
    v_t<ValueRow> m_values_foreach_rows_1;
    v_t<ValueRow> m_values_foreach_rows_2;

public:

    Smuvi(str_t name, const key_t& key, const value_t& value):
            m_inserter(name, InsertRow(key, value), {100, 2.0}),
            m_irank_world_to_iaccess_table(mpi::nrank(), ~0ul) {
        // index of the shared memory realm
        const auto ishmem = mpi::g_ishmems[mpi::irank()];
        m_access_tables.reserve(mpi::g_iranks_world_in_shmem_realms[ishmem].size());
        m_access_foreach_rows.reserve(m_access_tables.capacity());
        // iterate over the rank indices in this shmem realm and create an accessor table for each one
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            m_irank_world_to_iaccess_table[irank_world] = m_access_tables.size();
            m_access_tables.emplace_back(AccessRow(key), true);
            m_access_foreach_rows.emplace_back(m_access_tables.back().m_row);
        }
        /*
         * make sure enough elements are allocated for all the accessor tables
         */
        m_values_tables.reserve(m_access_tables.size());
        m_values_lookup_rows.reserve(m_access_tables.size());
        m_values_foreach_rows_1.reserve(m_access_tables.size());
        m_values_foreach_rows_2.reserve(m_access_tables.size());
        /*
         * loop over all the world-realm rank indices in this rank's shared memory region
         */
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            /*
             * create a new shared memory array of indices with sufficient elements to store all the entries associated
             * with all the keys sent to irank_world
             */
            m_values_tables.emplace_back(ValueRow(value), true);
            DEBUG_ASSERT_FALSE(m_values_tables.back().m_row.is_deref_valid(), "rows should be pointing to null");
            m_values_lookup_rows.emplace_back(m_values_tables.back().m_row);
            m_values_foreach_rows_1.emplace_back(m_values_tables.back().m_row);
            m_values_foreach_rows_2.emplace_back(m_values_tables.back().m_row);
        }
    }

    void insert(const key_t& key, const value_t& value) {
        // send a copy to one rank in each shared memory realm
        for (auto irank_dst: Distribution::one_irank_in_each_shmem_region(key)) {
            Table<InsertRow> &send = m_inserter.send(irank_dst);
            send.m_row.push_back_jump();
            send.m_row.m_key = key;
            send.m_row.m_value = value;
        }
    }

    struct AccessResult {
        const ValueRow& m_value_row;
        const uint_t m_index_end;

        operator bool () const {
            if (!m_value_row.is_deref_valid()) return false;
            return m_value_row.in_range(m_index_end);
        }

        bool operator==(const AccessResult& other) const {
            if (!m_value_row.is_deref_valid() && other.m_value_row.is_deref_valid()) return true;
            return m_value_row.index() == other.m_value_row.index() && m_index_end == other.m_index_end;
        }

        template<typename fn_t>
        void foreach(const fn_t& fn) const {
            functor::assert_prototype<void(const value_t&)>(fn);
            if (!*this) return;
            for (; m_value_row.in_range(m_index_end); ++m_value_row) fn(m_value_row.m_value);
        }
    };

    /**
     * lookup the value given a key
     */
    AccessResult access(const key_t& key, const ValueRow& value_iterator_row) const {
        // get the world rank index associated with the storage of this key
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        const auto itable = m_irank_world_to_iaccess_table[irank];
        const auto& accessor = m_access_tables[itable];
        const AccessRow lookup_row = accessor.lookup(key);
        if (!lookup_row) {
            // failed lookup
            value_iterator_row.select_null();
            return {value_iterator_row, ~0ul};
        }
        const uint_t ibegin = lookup_row.m_value_displ;
        const uint_t iend = ibegin + lookup_row.m_value_count;
        value_iterator_row.jump(ibegin);
        return {value_iterator_row, iend};
    }

    uint_t itable(const key_t& key) const {
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        return m_irank_world_to_iaccess_table[irank];
    }


    /**
     * lookup the value given a key
     */
    AccessResult access(const key_t& key) const {
        // get the world rank index associated with the storage of this key
        const auto itable = this->itable(key);
        const auto& value_row = m_values_lookup_rows[itable];
        return access(key, value_row);
    }

    template<typename fn_t>
    void foreach_value(const key_t& key, const fn_t& fn) const {
        functor::assert_prototype<void(const value_t&)>(fn);
        const auto itable = this->itable(key);
        const auto& value_row = m_values_foreach_rows_1[itable];
        access(key, value_row).foreach(fn);
    }

    /**
     * yield by call to fn_t each time a common value is found among the (ordered) values of the two keys
     */
    template<typename fn_t>
    void foreach_common_value(const key_t& key, const Smuvi<key_t, value_t>& other, const key_t& key_other, const fn_t& fn) const {
        functor::assert_prototype<void(const value_t&)>(fn);
        const auto itable_this = this->itable(key);
        const auto& value_row_this = m_values_foreach_rows_1[itable_this];
        const auto itable_other = other.itable(key_other);
        const auto& value_row_other = m_values_foreach_rows_2[itable_other];
        auto access_result_this = access(key, value_row_this);
        if (!access_result_this) return;
        auto access_result_other = access(key_other, value_row_other);
        if (!access_result_other) return;
        while (access_result_this && access_result_other) {
            if (value_row_this.m_value < value_row_other.m_value) ++value_row_this;
            else if (value_row_this.m_value > value_row_other.m_value) ++value_row_other;
            else fn(value_row_this.m_value);
        }
    }

    /**
     * overload in case both key-value sets are to be taken from the same SMUVI
     */
    template<typename fn_t>
    void foreach_common_value(const key_t& key_1, const key_t& key_2, const fn_t& fn) const {
        foreach_common_value(key_1, *this, key_2, fn);
    }

    template<typename fn_t>
    void foreach_value_pair(const key_t& key, const fn_t& fn, bool ordered = false, bool allow_eq = false) const {
        functor::assert_prototype<void(const value_t&, const value_t&)>(fn);
        const auto itable = this->itable(key);
        const auto& outer = m_values_foreach_rows_1[itable];
        const auto& inner = m_values_foreach_rows_2[itable];
        const AccessResult access_result = access(key, outer);
        if (!access_result) return;
        const auto ibegin = access_result.m_value_row.index();
        const auto iend = access_result.m_index_end;
        for (; outer.in_range(iend); ++outer) {
            if (ordered) {
                for (inner.jump(ibegin); inner.in_range(outer.index() + allow_eq); ++inner)
                    fn(outer.m_value, inner.m_value);
            }
            else {
                for (inner.jump(ibegin); inner.in_range(iend); ++inner) {
                    if (inner.index() != outer.index()) fn(outer.m_value, inner.m_value);
                }
            }
        }
    }

    template<typename fn_t>
    void foreach_key(const fn_t& fn) {
        functor::assert_prototype<void(const key_t&, AccessResult access_result)>(fn);
        DEBUG_ASSERT_EQ(m_access_tables.size(), m_values_tables.size(), "should have same number of access and values tables");
        for (auto itable = 0ul; itable < m_access_tables.size(); ++itable) {
            const auto& access_row = m_access_foreach_rows[itable];
            const auto& value_row = m_values_foreach_rows_1[itable];
            for (access_row.restart(); access_row; ++access_row) {
                const uint_t displ = access_row.m_value_displ;
                const uint_t count = access_row.m_value_count;
                value_row.jump(displ);
                fn(access_row.m_key, {value_row, displ + count});
            }
        }
    }

    void collate() {
        /*
         * first send the inserted key-index pairs to the receiving ranks
         */
        m_inserter.communicate();
        /*
         * a shared memory MappedTable is set up for each rank, so that they can be set up without dataraces by a single
         * rank, and thereafter accessed by all ranks in the shared memory region
         */
        MappedTable<AccessRow>& accessor = m_access_tables[m_irank_world_to_iaccess_table[mpi::irank()]];
        /*
         * initialize a private set of indices for each key sent to this rank. The sets are used to build up an
         * unordered list of the value element indices in the recv table associated with each key in the accessor,
         * which will later be copied into the shared memory m_entries arrays
         */
        const auto comp_recv_row_1 = m_inserter.recv().m_row;
        const auto comp_recv_row_2 = m_inserter.recv().m_row;
        auto comp_fn = [&](uint_t i, uint_t j) -> bool {
            comp_recv_row_1.jump(i);
            comp_recv_row_2.jump(j);
            return comp_recv_row_1.m_value < comp_recv_row_2.m_value;
        };
        v_t<std::set<uint_t, decltype(comp_fn)>> value_index_sets;

        /*
         * loop over the received rows
         */
        auto& recv_row = m_inserter.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            /*
             * lookup the received key in the accessor table, or insert it if this is the first instance
             */
            auto& accessor_row = accessor.lookup_or_insert(recv_row.m_key);
            /*
             * if there aren't enough value index vector sets for the current size of the accessor, allocate more
             */
            if (value_index_sets.size() < accessor.nrecord())
                value_index_sets.resize(accessor.nrecord(), std::set<uint_t, decltype(comp_fn)>(comp_fn));
            /*
             * get the indices vector to recv_row.m_key, and append the new index
             */
            value_index_sets[accessor_row.index()].insert(recv_row.index());
        }

        /*
         * a vector storing the number of entry_sets associated with each MPI rank
         */
        uintv_t nentries;
        {
            // the displacement from the beginning of the entry_sets
            uint_t displ = 0ul;
            /*
             * iterator initially pointing to the beginning of the entry sets (i.e. the first key in the accessor
             * created by this rank)
             */
            auto value_index_set_it = value_index_sets.cbegin();
            // loop over all the rows in the accessor created by this rank
            auto& accessor_row = accessor.m_row;
            for (accessor_row.restart(); accessor_row; ++accessor_row) {
                /*
                 * set the displacement that will denote the index in the entry_sets array at which the entry_sets
                 * corresponding to the current row in the accessor (i.e. the key) begin
                 */
                accessor_row.m_value_displ = displ;
                /*
                 * set the number of entry_sets corresponding to the key
                 */
                accessor_row.m_value_count = value_index_set_it->size();
                /*
                 * increment the rank-private displ count
                 */
                displ += accessor_row.m_value_count;
                /*
                 * advance the entry iterator
                 */
                ++value_index_set_it;
            }
            DEBUG_ASSERT_TRUE(value_index_set_it == value_index_sets.cend(), "should have iterated through all entry sets");
            /*
             * gather the final displs (i.e. the total number of unique values across all keys sent to this rank)
             */
            mpi::all_gather(displ, nentries);
        }

        /*
         * clear the entries shared memory arrays in case this is not the first call to collate
         */
        for (auto& values_table: m_values_tables) values_table.clear();


        /*
         * now write the sets of entries stored privately on this rank to the shared memory arrays accessible by all
         * ranks on the same shared memory region
         */
        {
            /*
             * get a reference to the values array updated by this rank
             */
            auto& values_table = m_values_tables[m_irank_world_to_iaccess_table[mpi::irank()]];
            DEBUG_ASSERT_TRUE(values_table.empty(), "value tables should be been cleared");
            values_table.push_back(nentries[mpi::irank()]);
            uint_t ielement = 0;
            auto& value_row = values_table.m_row;
            for (const auto& value_index_set : value_index_sets) {
                /*
                 * loop over all value indices in all sets, copying the indexed values contiguously to the values table
                 */
                for (const auto& i: value_index_set) {
                    recv_row.jump(i);
                    value_row.jump(ielement++);
                    value_row.m_value = recv_row.m_value;
                }
            }
        }
        /*
         * now the accessor tables are updated with the keys, and the entries arrays are updated with their
         * corresponding sorted and unique index values
         */
    }
};

#endif //M7_SMUVI_H
