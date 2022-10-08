//
// Created by Robert J. Anderson on 08/05/2021.
//

#ifndef M7_TABLEBASE_H
#define M7_TABLEBASE_H

#include <set>
#include <M7_lib/io/Logging.h>
#include <numeric>
#include "Buffer.h"

struct RowTransfer {
    /**
     * make rows to be sent contiguous in memory
     */
    Buffer m_send_buffer, m_recv_buffer;
    Buffer::Window m_send_bw, m_recv_bw;
    /**
     * need to get a unique pair of P2P tag from global variables
     */
    const int m_nrec_p2p_tag = mpi::new_p2p_tag();
    const int m_irecs_p2p_tag = mpi::new_p2p_tag();

    RowTransfer(str_t name) :
            m_send_buffer("Outward transfer buffer", 1),
            m_recv_buffer("Inward transfer buffer", 1) {
        logging::info("Initializing row send/recv buffers for table \"{}\"", name);
        logging::debug("P2P tag for number of row indices to transfer for \"{}\": {}", name, m_nrec_p2p_tag);
        logging::debug("P2P tag for array of row indices to transfer for \"{}\": {}", name, m_irecs_p2p_tag);
        m_send_buffer.append_window(&m_send_bw);
        m_recv_buffer.append_window(&m_recv_bw);
    }
};

struct RowProtector;

/**
 * Base class for all row-contiguous data in the program.
 *
 * The data layout of a Table is defined by a class derived from Row, but this templated dependence is not included
 * here in this base class. While the fine detail of field offsets etc is not needed here, we need to know the basics,
 * i.e. what amount of memory does a single row require. This is always padded to a whole number of system words.
 * A system word is sizeof(size_t) bytes in length
 *
 * So as to avoid confusion between the Row class (and its subclasses) and the data buffer slice used for writing
 * and reading the fields of a Row - the latter are termed "slots". Slots can be filled with a "record", that is the
 * raw data contents managed by a Row-derived object, or they can be empty, in which case they are denoted as "free
 * slots".
 *
 * The contents of the Table are ultimately stored in a Buffer, which is just a wrapper for a dynamically
 * allocated array of buf_t elements. Buffers can be shared by many Tables though, so a table is instead
 * given access to a "window" of that Buffer. All tables which have windows on the same buffer have the ability
 * to resize, and in response to this, the Buffer will ensure that table records are moved to new positions in the
 * resized buffer, and that each of the other windows are pointed to the beginning of that data.
 *
 * The number of slots is determined by the size of the currently-allocated buffer window. This is analogous to the
 * capacity of a std::vector. The "high water mark" (m_hwm) is the analog of the size() of a std::vector. The number of
 * records is given by the difference of the high water mark and the number of free slots
 *
 * the m_slot_size member of m_bw is the definitive length of the Table's row in bytes, and the Row class ensures
 * that this is always an integer multiple of the system word length
 */
struct TableBase {
    /**
     * this Table's portion of the raw data buffer possibly shared by many Tables
     */
    Buffer::Window m_bw;
    /**
     * "high water mark" is the number of slots in use (free slots and records)
     */
    uint_t m_hwm = 0ul;
    /**
     * indices of freed slots, i.e. those slots below the high water mark that have been cleared
     */
    std::stack<uint_t> m_freed_slots;
    /**
     * a constant-time lookup for whether a given slot is freed i.e. true if it appears in m_freed_slots
     */
    std::vector<bool> m_is_freed_slot;


    /**
     * buffered space for communication in the event of rank reallocation, instantiate on first transfer if required
     */
    std::unique_ptr<RowTransfer> m_transfer = nullptr;
    /**
     * if a record is to be protected, its slot index must appear as a key in the following map.
     * the value associated with the key is the "protection level" of the record. if one object requires that a record i
     * not be deleted, then {i: 1} appears in the map. If a second object requires that record i not be deleted,
     * then the key-value pair is looked up, then the value is incremented so that {i: 2} appears in the map. If then
     * the first object releases its protection of record i, the key-value pair then returns to {i: 1}, so that the
     * record will still not be deleted even though the first dependent object has no objection. The record will remain
     * protected (indelible) until the second object releases, and the corresponding key-value pair record is deleted
     * from the map
     */
    mutable std::map<uint_t, uint_t> m_protected_records;

    TableBase(uint_t slot_size);

protected:
    /**
     * the copy assignment is not exposed as public, since otherwise the TableBase could be set equal to an incompatible
     * TableBase. Better and safer that the row-templated subclasses implement the public interface
     * @param other
     *  value of TableBase to which this one will be set
     * @return
     *  ref to this
     */
    TableBase& operator=(const TableBase& other) {
        if (other.nslot()) {
            resize(other.nslot(), 0.0);
            m_bw = other.m_bw;
        }
        m_hwm = other.m_hwm;
        m_freed_slots = other.m_freed_slots;
        m_protected_records = other.m_protected_records;
        return *this;
    }

public:
    TableBase(const TableBase& other);

    /**
     * the total number of slots in the BufferWindow. Think of this as the Table's capacity by analogy to std::vector
     */
    uint_t nslot() const {
        return m_bw.m_nslot;
    }

    /**
     * strict equality: other table may have the same contents but a different buffer size, in this case the tables are
     * not considered equal. main use for this method is to check for successful copy operations
     * @param other
     *  table to compare against
     * @return
     *  true if the other table is identical in value to this one
     */
    bool operator==(const TableBase& other) const {
        if (this==&other) return true;
        if (nslot() != other.nslot()) return false;
        if (m_hwm!=other.m_hwm) return false;
        return std::memcmp(begin(), other.begin(), m_hwm * slot_size()) == 0;
    }

    bool operator!=(const TableBase& other) const {
        return !(*this==other);
    }
    /**
     * the size requirement of a single row in bytes (always an integer number of system words)
     */
    uint_t slot_size() const {
        return m_bw.m_slot_size;
    }

    /**
     * @return
     *  pointer to the first data word of the BufferWindow
     */
    buf_t *begin();

    /**
     * @return
     * const pointer to the first data word of the BufferWindow
     */
    const buf_t *begin() const;

    /**
     * @param islot
     *  slot index
     * @return
     *  pointer to the beginning of the indexed slot
     */
    buf_t *begin(uint_t islot);

    /**
     * @param islot
     *  slot index
     * @return
     *  const pointer to the beginning of the indexed slot
     */
    const buf_t *begin(uint_t islot) const;

    /**
     * Associate the table with a buffer by assigning the table an available BufferWindow
     * @param buffer
     *  pointer to buffer
     */
    void set_buffer(Buffer *buffer);

    /**
     * increase the high water mark, expand first if necessary
     * @param nslot
     *  number of new slots to create
     * @return
     *  index of the newly-accessible record
     */
    uint_t push_back(uint_t nslot=1ul);

    /**
     * If there are indices on the m_free_slots stack: pop one and use it, else: push_back
     * @return
     *  index of free slot
     */
    uint_t get_free_slot();

    /**
     * clear the entire table only if it contains no protected records, else fatal error
     */
    virtual void clear();

    /**
     * clear the byte string of the indexed slot and designate it "freed". If it is already freed or protected, a fatal
     * error is thrown
     * @param islot
     *  index of slot to free
     */
    virtual void free(uint_t islot);

    bool empty() const;

    bool is_freed(uint_t islot) const;

    /**
     * @return
     *  size of the buffer window in bytes
     */
    uint_t bw_size() const;

    /**
     * call the resize method on the buffer window and reflect the reallocation in m_nslot
     * @param nrec
     *  minimum number of records in the new buffer.
     */
    virtual void resize(uint_t nrec, double factor=-1.0);

    /**
     * resize based on the number of additional records required beyond those currently allocated
     * @param nrec
     *  minimum number of new records.
     */
    void expand(uint_t nrec, double factor=-1.0);

    /**
     * simply call clear on each indexed record
     * @param irecs
     *  all record indices marked for erasure
     */
    void free_records(const uintv_t& irecs);

    /**
     * in some derived classes, there is more to adding new records than simply copying their contents into the hwm. In
     * these cases we need to call a method after the copy in order to maintain the integrity of those data structures
     * see MappedTable.h for more details
     * @param iinsert
     *  record index of the already-inserted data
     */
    virtual void post_insert(uint_t /*iinsert*/){}

    /**
     * If we have copied a block of records contiguously into a table with non-trivial post-insert obligations, we need
     * to call post_insert for each copied record
     * @param ibegin
     *  record index of beginning of the already-inserted data
     * @param iend
     *  record index of end of the already-inserted data
     */
    void post_insert_range(uint_t ibegin = 0ul, uint_t iend = ~0ul) {
        if (iend == ~0ul) iend = m_hwm;
        for (uint_t i = ibegin; i < iend; ++i) post_insert(i);
    }

    /**
     * function pointer type for the callback associated with row transfers.
     * see RankAllocator.h
     */
    typedef std::function<void(const uintv_t& , uint_t, uint_t)> transfer_cb_t;
    /**
     * function pointer type for the callback associated with receipt of a single row in a transfer operation
     * see RankAllocator.h
     */
    typedef std::function<void(uint_t)> recv_cb_t;

    /**
     * insert all records held in the buffer window into this Table, then call all callbacks if any.
     * @param recv
     *  received buffer containing transferred records contiguously, and densely - no zeroed records.
     * @param nrec
     *  number of records being inserted
     * @param callbacks
     *  functions to call each time a received record is processed
     */
    virtual void insert_records(const Buffer::Window& recv, uint_t nrec, const std::list<recv_cb_t>& callbacks);

    /**
     * By P2P MPI communication, send the records identified in the first arg from irank_send to irank_recv.
     * Called on all ranks, but the first arg is only respected on the irank_send rank.
     * @param irecs
     *  record indices on irank_send to be found, copied to contiguous buffer and transferred
     * @param irank_send
     *  MPI rank index from which the records are being transferred
     * @param irank_recv
     *  MPI rank index to which the records are being transferred
     * @param callbacks
     *  functions to call on irank_recv each time a received record is processed
     */
    void transfer_records(const uintv_t& irecs, uint_t irank_send, uint_t irank_recv,
                          const std::list<recv_cb_t>& callbacks = {});

    /**
     * copy a single record from another "source" table
     * @param src
     *  source table which must have the same record size. should also have the same overall data layout but this is not
     *  verified at this level
     * @param islot_src
     *  record in the source table to be copied
     * @param islot_dst
     *  record in this table to which the source record is copied
     */
    void copy_record_in(const TableBase& src, uint_t islot_src, uint_t islot_dst);

    /**
     * swap record contents via std::swap_ranges
     * @param irec
     *  index of record to be swapped
     * @param jrec
     *  index of other record to be swapped
     */
    void swap_records(uint_t irec, uint_t jrec);

    /**
     * increase protection level of the record. insert key-value pair first if non-existent
     * @param islot
     *  slot index of the record to protect
     */
    void protect(uint_t islot) const {
        DEBUG_ASSERT_LT(islot, m_hwm, "slot index OOB");
        DEBUG_ASSERT_FALSE(m_is_freed_slot[islot], "cannot protect a freed slot");
        auto it = m_protected_records.find(islot);
        if (it == m_protected_records.end()) m_protected_records.insert({islot, 1ul});
        else ++it->second;
    }

    /**
     * @param islot
     *  slot index
     * @return
     *  number of protect(islot) calls - number of release(islot) calls
     */
    uint_t protection_level(uint_t islot) const {
        DEBUG_ASSERT_LT(islot, m_hwm, "slot index OOB");
        auto it = m_protected_records.find(islot);
        if (it == m_protected_records.end()) return 0ul;
        else return it->second;
    }

    /**
     * @param islot
     *  slot index of record
     * @return
     *  true if the protection level is non-zero
     */
    bool is_protected(uint_t islot) const {
        return protection_level(islot);
    }

    /**
     * remove protection from the indexed record
     * @param islot
     *  slot index of protected record
     */
    void release(uint_t islot) const {
        DEBUG_ASSERT_TRUE(protection_level(islot), "can't release an unprotected row");
        auto it = m_protected_records.find(islot);
        if (it->second==1ul) m_protected_records.erase(it);
        --it->second;
    }

    /**
     * debugging method to ensure consistency of the two structures keeping track of the free status of slots
     * @return
     *  true is m_free_slots member is consistent with m_is_free_slot
     */
    bool freed_slots_consistent() const;

    /**
     * "Location" class which describes the location of a record in a distributed table i.e. by a local slot index and
     * an MPI rank index
     */
    struct Loc {
        /**
         * rank and slot indices
         */
        const uint_t m_irank, m_islot;

        Loc(uint_t irank, uint_t islot);

        /**
         * @return
         *  true if the location is anywhere i.e. has a valid rank index
         */
        operator bool() const;

        operator uintv_t() const {
            if (is_mine()) return {m_islot};
            return {};
        }

        /**
         * @return
         *  true if MPI rank index matches bcast-shared rank of identified slot
         */
        bool is_mine() const;

        bool operator==(const Loc& other);

        bool operator!=(const Loc& other);
    };

    /**
     * When Field-based data structure is introduced in the derived classes, this method is capable of displaying
     * human-readable columns. here though, we don't know the data layout, so return slot bytes as integers
     * @param ordering
     *  optional indices to reorder table on the fly as it is printed, without the need to physically reorder slots
     * @return
     *  string representing table's contents
     */
    virtual str_t to_string(const uintv_t* ordering = nullptr) const;

    /**
     * gather contents of another table over all MPI ranks into this table on all ranks
     * @param src
     *  table whose slots are to be gathered.
     */
    virtual void all_gatherv(const TableBase& src);

    /**
     * gather contents of another table over all MPI ranks into this table on the given rank
     * @param src
     *  table whose slots are to be gathered.
     * @param irank
     *  index of the only rank in the communicator to receive the data from src
     */
    virtual void gatherv(const TableBase& src, uint_t irank = 0ul);

    /**
     * @return
     *  true if any of the associated RowProtectors are protecting any records of this table
     */
    bool is_protected() const;

    /**
     * @return
     *  number of slots below the high water mark that aren't free
     */
    uint_t nrecord() const;

    str_t name() const {
        return m_bw.name();
    }

    void rename(str_t name) const {
        m_bw.rename(name);
    }
};


#endif //M7_TABLEBASE_H
