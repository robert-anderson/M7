//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_SMUVI_H
#define M7_SMUVI_H

#include "MappedTable.h"
#include "BufferedTable.h"
#include "M7_lib/communication/SendRecv.h"
#include "M7_lib/communication/Distribution.h"


struct SmuviEntries {
    /**
     * pointer to first element mapped to by a given key
     */
    const size_t* m_cbegin;
    /**
     * pointer to first element after m_cbegin not mapped to by the given key
     */
    const size_t* m_cend;

    operator bool () const {
        return m_cbegin && m_cend && (m_cbegin != m_cend);
    }

    template<typename fn_t>
    void foreach(const fn_t& fn) const {
        functor::assert_prototype<void(uint_t)>(fn);
        for (auto ptr = m_cbegin; ptr != m_cend; ++ptr) fn(*ptr);
    }
};

struct SmuviEntriesWithIndex {
    const SmuviEntries m_entries;
    const size_t m_key_index;
    operator bool() const {
        return m_entries;
    }
};

/**
 * Shared Memory Unordered map to Vectors of Indices(/Items) "SMUVI"
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


template<typename key_t>
class Smuvi {

    /**
     * Row type for the inserter tables
     */
    struct InsertRow : Row {
        key_t m_key;
        field::Number<uint_t> m_entry;

        template<typename ...Args>
        InsertRow(const Args&... key_ctor_args): m_key(this, key_ctor_args...), m_entry(this){}

        key_t &key_field() {
            return m_key;
        };
    };
    /**
     * Row type for the accessor tables
     */
    struct AccessRow : Row {
        key_t m_key;
        field::Number<uint_t> m_entry_count;
        field::Number<uint_t> m_entry_displ;

        template<typename ...Args>
        AccessRow(const Args&... key_ctor_args):
                m_key(this, key_ctor_args...), m_entry_count(this), m_entry_displ(this){}

        key_t &key_field() {
            return m_key;
        };
    };
    /**
     * send and receive tables for the key-value pair insertion
     */
    send_recv::BasicSend<InsertRow> m_inserter;
    uintv_t m_irank_world_to_iaccessor;
    /**
     * one shared-memory MappedTable for each rank in the shared memory region
     */
    v_t<buffered::MappedTable<AccessRow>> m_accessors;
    /**
     * a single flat storage for all the entries associated with all keys in the shared memory region
     */
    v_t<SharedArray<uint_t>> m_entries;

public:

    template<typename ...Args>
    Smuvi(str_t name, const Args&... key_ctor_args):
            m_inserter(name, InsertRow(key_ctor_args...), {100, 2.0}),
            m_irank_world_to_iaccessor(mpi::nrank(), ~0ul) {
        // index of the shared memory realm
        const auto ishmem = mpi::g_ishmems[mpi::irank()];
        // iterate over the rank indices in this shmem realm and create an accessor table for each one
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            m_irank_world_to_iaccessor[irank_world] = m_accessors.size();
            m_accessors.emplace_back(AccessRow(key_ctor_args...), true);
        }
    }

    void insert(const key_t& key, const uint_t& entry) {
        // send a copy to one rank in each shared memory realm
        for (auto irank_dst: Distribution::one_irank_in_each_shmem_region(key)) {
            Table<InsertRow> &send = m_inserter.send(irank_dst);
            send.m_row.push_back_jump();
            send.m_row.m_key = key;
            send.m_row.m_entry = entry;
        }
    }

    /**
     * lookup entries (a vector of indices) and the key index given a key
     */
    SmuviEntriesWithIndex access_by_key(const key_t& key) const {
        // get the world rank index associated with the storage of this key
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        auto& accessor = m_accessors[m_irank_world_to_iaccessor[irank]];
        auto& entries = m_entries[m_irank_world_to_iaccessor[irank]];
        const AccessRow& lookup_row = accessor.lookup(key);
        if (!lookup_row) return {{nullptr, nullptr}, ~0ul}; // failed lookup
        const auto cbegin = entries.cbegin() + lookup_row.m_entry_displ;
        // todo: rja, make index global wrt shared memory region offsets
        return {{cbegin, cbegin + lookup_row.m_entry_count}, lookup_row.index()};
    }

    struct SmuviEntriesWithKey {
        const SmuviEntries m_entries;
        const AccessRow& m_key_row;
        operator bool() const {
            return m_entries;
        }
    };
    /**
     * lookup entries and the corresponding key given an index
     */
    SmuviEntriesWithKey access_by_index(size_t index) const {
        // todo: rja, make shmem-region aware
        // get the world rank index associated with the storage of this key
        auto irank = 0ul;
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        auto& accessor = m_accessors[m_irank_world_to_iaccessor[irank]];
        auto& entries = m_entries[m_irank_world_to_iaccessor[irank]];
        const auto& lookup_row = accessor.m_lookup_row;
        if (index >= accessor.nrow_in_use()) return {{nullptr, nullptr}, lookup_row}; // failed lookup, out of bounds
        lookup_row.jump(index);
        const auto cbegin = entries.cbegin() + lookup_row.m_entry_displ;
        // todo: rja, make index global wrt shared memory region offsets
        return {{cbegin, cbegin + lookup_row.m_entry_count}, lookup_row};
    }

    template<typename fn_t>
    void foreach_entry_by_key(const key_t& key, const fn_t& fn) const {
        access_by_key(key).m_entries.foreach(fn);
    }

    template<typename fn_t>
    void foreach_entry_by_index(uint_t index, const fn_t& fn) const {
        access_by_index(index).m_entries.foreach(fn);
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
        MappedTable<AccessRow>& accessor = m_accessors[m_irank_world_to_iaccessor[mpi::irank()]];
        /*
         * initialize a private vector of sets, once set for each key sent to this rank. The sets are used to build
         * up an ordered list of the unique entry_sets associated with each key in the accessor, which will later be copied
         * into the shared memory m_entries arrays
         */
        v_t<std::set<uint_t>> entry_sets;
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
             * if there aren't enough entry sets for the current size of the accessor, allocate more
             */
            if (entry_sets.size() < accessor.nrecord()) entry_sets.resize(accessor.nrecord());
            /*
             * get the set of entry_sets corresponding to recv_row.m_key, and insert the new entry if it doesn't exist
             */
            entry_sets[accessor_row.index()].insert(recv_row.m_entry);
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
            auto entry_sets_it = entry_sets.cbegin();
            // loop over all the rows in the accessor created by this rank
            auto& accessor_row = accessor.m_row;
            for (accessor_row.restart(); accessor_row; ++accessor_row) {
                /*
                 * set the displacement that will denote the index in the entry_sets array at which the entry_sets
                 * corresponding to the current row in the accessor (i.e. the key) begin
                 */
                accessor_row.m_entry_displ = displ;
                /*
                 * set the number of entry_sets corresponding to the key
                 */
                accessor_row.m_entry_count = entry_sets_it->size();
                /*
                 * increment the rank-private displ count
                 */
                displ += accessor_row.m_entry_count;
                /*
                 * advance the entry iterator
                 */
                ++entry_sets_it;
            }
            DEBUG_ASSERT_TRUE(entry_sets_it == entry_sets.end(), "should have iterated through all entry sets");
            /*
             * gather the final displs (i.e. the total number of entries across all keys sent to this rank)
             */
            mpi::all_gather(displ, nentries);
        }

        /*
         * clear the entries shared memory arrays in case this is not the first call to collate
         */
        m_entries.clear();
        /*
         * make sure enough elements are allocated for all the accessor tables
         */
        m_entries.reserve(m_accessors.size());
        /*
         * get the shared memory region index of this rank
         */
        const auto ishmem = mpi::g_ishmems[mpi::irank()];
        /*
         * loop over all the world-realm rank indices in this rank's shared memory region
         */
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            /*
             * create a new shared memory array of indices with sufficient elements to store all the entries associated
             * with all the keys sent to irank_world
             */
            m_entries.emplace_back(SharedArray<uint_t>(nentries[irank_world]));
        }

        /*
         * now write the sets of entries stored privately on this rank to the shared memory arrays accessible by all
         * ranks on the same shared memory region
         */
        {
            /*
             * get a reference to the entries array updated by this rank
             */
            auto& entries_array = m_entries[m_irank_world_to_iaccessor[mpi::irank()]];
            uint_t ientry = 0;
            for (const auto& entry_set : entry_sets) {
                /*
                 * loop over all values in all entry sets, copying them contiguously to the entries array
                 */
                for (const auto& entry: entry_set) entries_array.set_(ientry++, entry);
            }
        }
        /*
         * now the accessor tables are updated with the keys, and the entries arrays are updated with their
         * corresponding sorted and unique index values
         */
    }
};

#endif //M7_SMUVI_H
