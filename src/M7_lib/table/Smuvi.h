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
 * Shared Memory Unordered map to Vectors of Indices(/Items) "SMUVI"
 * The SMUVI is a type of parallel hash map from a key domain to a variable-length vector of indices which has fewer
 * requirements than a parallelised generalization of a map from type T to a arbitrarily sized array of type U elements,
 * expressed using STL containers as std::unordered_map<T, std::vector<unsigned long>>.
 *
 * Primary among which is that the map is never simultaneously modified and read; it is at any given time in one of two
 * "insert" or "access"
 *
 * The required data structure need not support some operations commonly implemented in serial and concurrent hash maps:
 *  - Reassignment of existing elements. Newly inserted scalar values are appended to the existing list of values
 *    corresponding to the key
 *  - Deletion. No key or element of a value array may be removed after it has been added
 *  - Dynamic resizing of the "access" key array.
 *
 * The SMUVI behaves as a hash table which can be inserted to in parallel, then __collated__ with minimal
 * synchronisation overhead, to then be accessed for reading in constant time among all ranks in a communicator with a
 * common shared memory window.
 *
 * In the insert state, each rank initialises a private hash map and this is filled according to the workload of
 * insertions allotted to the rank. For access, shared-memory implementations of hash maps usually need to allow for the
 * possibility that a bin of keys can be written to or deleted from while the read operation is being carried out, and
 * as such a lock must be acquired before a read can take place, which results in sychronisation overhead. However, if
 * the SMUVI is in access state, each rank in a shared memory realm can lookup elements without the need to lock a
 * region of the buffer.
 *
 * The SMUVI state is transformed from insert to access in the collate operation. This proceeds by collecting the keys
 * from each insertion map into a local array which is then sorted. These sorted local keys are then gathered
 * contiguously into a shared memory buffer.
 */


template<typename key_t>
struct Smuvi {

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

private:
    const MappedTable<AccessRow>& accessor(const key_t& key) const {
        // get the world rank index associated with the storage of this key
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        return m_accessors[m_irank_world_to_iaccessor[irank]];
    }

    const SharedArray<uint_t>& entries(const key_t& key) const {
        // get the world rank index associated with the storage of this key
        auto irank = Distribution::irank_in_shmem_region(key);
        DEBUG_ASSERT_LT(irank, ~0ul, "MPI rank should be assigned an allocated accessor");
        return m_entries[m_irank_world_to_iaccessor[irank]];
    }

public:
    void insert(const key_t& key, const uint_t& entry) {
        // send a copy to one rank in each shared memory realm
        for (auto irank_dst: Distribution::one_irank_in_each_shmem_region(key)) {
            Table<InsertRow> &send = m_inserter.send(irank_dst);
            send.m_row.push_back_jump();
            send.m_row.m_key = key;
            send.m_row.m_entry = entry;
        }
    }

    uint_t nitem(const key_t& key) const {
        const MappedTable<AccessRow>& accessor = this->accessor(key);
        auto& entries_vec = m_entries[m_irank_world_to_iaccessor[mpi::irank()]];
        const AccessRow& lookup_row = accessor.lookup(key);
        if (!lookup_row) return ~0ul; // failed lookup
        return lookup_row.m_entry_count;
    }

    const uint_t* item_cbegin(const key_t& key) const {
        const MappedTable<AccessRow>& accessor = this->accessor(key);
        const AccessRow& lookup_row = accessor.lookup(key);
        if (!lookup_row) return nullptr; // failed lookup
        return entries(key).cbegin() + lookup_row.m_entry_displ;
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


        m_entries.clear();
        m_entries.reserve(m_accessors.size());

        const auto ishmem = mpi::g_ishmems[mpi::irank()];
        for (auto& irank_world: mpi::g_iranks_world_in_shmem_realms[ishmem]) {
            m_entries.emplace_back(SharedArray<uint_t>(nentries[irank_world]));
        }


        {
            auto& accessor_row = accessor.m_row;
            auto& entries_vec = m_entries[m_irank_world_to_iaccessor[mpi::irank()]];
            uint_t iitem = 0;
            for (auto items : entries) {
                for (auto item: items) entries_vec.set_(iitem++, item);
            }
        }
    }
};

#endif //M7_SMUVI_H
