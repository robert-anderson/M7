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

    struct InsertRow : Row {
        key_t m_key;
        field::Number<uint_t> m_entry;

        template<typename ...Args>
        InsertRow(const Args&... key_ctor_args): m_key(this, key_ctor_args...), m_entry(this){}

        key_t &key_field() {
            return m_key;
        };
    };

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

    send_recv::BasicSend<InsertRow> m_inserter;
    uintv_t m_irank_world_to_iaccessor;
    v_t<buffered::MappedTable<AccessRow>> m_accessors;
    v_t<SharedArray<uint_t>> m_entries;
    v_t<uint_t> m_rank_offsets;
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

    SmuviEntries access(const key_t& key) const {
        // get the world rank index associated with the storage of this key
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        auto& accessor = m_accessors[m_irank_world_to_iaccessor[irank]];
        auto& entries = m_entries[m_irank_world_to_iaccessor[irank]];
        const AccessRow& lookup_row = accessor.lookup(key);
        if (!lookup_row) return {nullptr, nullptr}; // failed lookup
        const auto cbegin = entries.cbegin() + lookup_row.m_entry_displ;
        return {cbegin, cbegin + lookup_row.m_entry_count};
    }

    template<typename fn_t>
    void foreach_entry(const key_t& key, const fn_t& fn) const {
        access(key).foreach(fn);
    }

    void collate() {
        m_inserter.communicate();
        MappedTable<AccessRow>& accessor = m_accessors[m_irank_world_to_iaccessor[mpi::irank()]];
        v_t<std::set<uint_t>> entries;

        auto& recv_row = m_inserter.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            auto& accessor_row = accessor.lookup_or_insert(recv_row.m_key);
            if (entries.size() < accessor.nrecord()) entries.resize(accessor.nrecord());
            entries[accessor_row.index()].insert(recv_row.m_entry);
        }

        uintv_t nentries;
        {
            uint_t displ = 0ul;
            auto& accessor_row = accessor.m_row;
            auto entry_it = entries.cbegin();
            for (accessor_row.restart(); accessor_row; ++accessor_row) {
                accessor_row.m_entry_displ = displ;
                accessor_row.m_entry_count = entry_it->size();
                displ += accessor_row.m_entry_count;
                ++entry_it;
            }
            mpi::all_gather(displ, nentries);
        }

        m_rank_offsets.clear();
        m_rank_offsets.push_back(0);
        for (auto it = m_accessors.cbegin(); it != m_accessors.cend()-1; ++it) {
            m_rank_offsets.push_back(m_rank_offsets.back() + it->nrow_in_use());
        }


        m_entries.clear();
        m_entries.reserve(m_accessors.size());

        const auto ishmem = mpi::g_ishmems[mpi::irank()];
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            m_entries.emplace_back(SharedArray<uint_t>(nentries[irank_world]));
        }

        {
            auto& entries_vec = m_entries[m_irank_world_to_iaccessor[mpi::irank()]];
            uint_t iitem = 0;
            for (auto& items : entries) {
                for (auto item: items) entries_vec.set_(iitem++, item);
            }
        }


    }
};

#endif //M7_SMUVI_H
