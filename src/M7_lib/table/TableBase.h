//
// Created by Robert J. Anderson on 08/05/2021.
//

#ifndef M7_TABLEBASE_H
#define M7_TABLEBASE_H

#include <set>
#include <M7_lib/io/Logging.h>
#include <numeric>
#include "Buffer.h"

/**
 * this forward declaration is needed so that SendRecv can be declared a friend and can therefore modify the size of the
 * usable portion of a Table
 */
template<typename row_t, typename send_table_t> class SendRecv;

/**
 * Base class for all row-contiguous data in the program.
 *
 * The data layout of a Table is defined by a class derived from Row, but this templated dependence is not included
 * here in this base class. While the fine detail of field offsets etc is not needed here, we need to know the basics,
 * i.e. what amount of memory does a single row require. This is always padded to a whole number of system words.
 * A system word is sizeof(size_t) bytes in length.
 *
 * Here some terminology is explained:
 *
 *  record
 *      a slice of the Table's BufferWindow associated with the storage of a Row
 *  cleared
 *      a row is cleared if its record is equal to a null buffer of the same length
 *  capacity
 *      the number of rows currently allocated
 *  high-water mark (HWM)
 *      the number of rows currently "in use", i.e. those that may be read and modified by a Row object.
 *  freed
 *      freeing a record is to clear it and mark it for reuse so that new entries can make use of the records left empty
 *      by clearing. The Row-wise clearing method is not made public, since all row deletions should be freeing so as
 *      not to leak recyclable records
 *  push back
 *      the act of increasing the HWM by a number (usually 1) of cleared (but not freed) records
 *  resize / expand
 *      modifies the size of the buffer window, and thereby the capacity. such calls reallocate the buffer and copy
 *      records below the HWM. Programmers should beware that such calls will invalidate the m_begin member of any rows
 *      associated with the resized Table.
 *  protected
 *      a protected record may not be removed by clearing
 *  protection level
 *      the number of reasons for which a record is protected (e.g. reference and semistochastic constituent =>
 *      protection level of 2)
 *  unprotect
 *      decrement the protection level of a record
 *
 *
 * The contents of the Table are ultimately stored in a Buffer, which is just a wrapper for a dynamically
 * allocated array of buf_t elements. Buffers can be shared by many Tables though, so a table is instead
 * given access to a "window" of that Buffer. All tables which have windows on the same buffer have the ability
 * to resize, and in response to this, the Buffer will ensure that table records are moved to new positions in the
 * resized buffer, and that each of the other windows are pointed to the beginning of that data.
 *
 * The number of records is determined by the size of the currently-allocated buffer window. This is analogous to the
 * capacity of a std::vector. The "nrow_in_use" method is the analog of the size() of a std::vector.
 *
 * the m_row_size member of m_bw is the definitive length of the Table's row in bytes, and the Row class ensures that
 * this is always an integer multiple of the system word length
 */
struct TableBase {

    template<typename row_t, typename send_table_t>
    friend class SendRecv;

protected:
    /**
     * this Table's portion of the raw data buffer possibly shared by many Tables
     */
    Buffer::Window m_bw;
    /**
     * indices of freed rows
     */
    std::stack<uint_t> m_freed_rows;
    /**
     * a constant-time lookup for whether a given row is freed i.e. true if it appears in m_freed_rows
     */
    std::vector<bool> m_is_freed_row;
    /**
     * retain a buffer of zeros for fast memcmp to determine whether a row is clear
     */
    const v_t<buf_t> m_null_row_string;

    /**
     * buffered space for communication in the event of rank reallocation, instantiate on first transfer if required
     */
//    std::unique_ptr<RowTransfer> m_transfer = nullptr;
    /**
     * if a record is to be protected, its index must appear as a key in the following map. the value associated with
     * the key is the "protection level" of the record. if one object requires that a record i not be deleted, then
     * {i: 1} appears in the map. If a second object requires that record i not be deleted, then the key-value pair is
     * looked up, then the value is incremented so that {i: 2} appears in the map. If then the first object releases its
     * protection of record i, the key-value pair then returns to {i: 1}, so that the record will still not be deleted
     * even though the first dependent object has no objection. The record will remain protected (indelible) until the
     * second object releases, and the corresponding key-value pair record is deleted from the map
     */
    mutable std::map<uint_t, uint_t> m_protected_rows;

    uint_t row_index(const buf_t* row_ptr) const {
        DEBUG_ASSERT_TRUE(row_ptr, "row pointer is undefined");
        const auto nbyte = std::distance<const buf_t*>(m_bw.cbegin(), row_ptr);
        DEBUG_ASSERT_GE(nbyte, 0l, "pointer is before the table beginning");
        DEBUG_ASSERT_FALSE(nbyte % m_bw.m_row_size, "pointer does not point to beginning of a row");
        return uint_t(nbyte / m_bw.m_row_size);
    }

    /**
     * the copy assignment is not exposed as public, since otherwise the TableBase could be set equal to an incompatible
     * TableBase. Better and safer that the row-templated subclasses implement the public interface
     * @param other
     *  value of TableBase to which this one will be set
     * @return
     *  ref to this
     */
    TableBase& operator=(const TableBase& other) {
        if (other.capacity()) {
            resize(other.capacity(), 0.0);
            m_bw = other.m_bw;
        }
        m_freed_rows = other.m_freed_rows;
        m_protected_rows = other.m_protected_rows;
        return *this;
    }

    /**
     * clear the indexed row only if it is not protected, else fatal error. This method is not public because when
     * erasure of a record is required, it should be marked as "free", not merely set to a null buffer
     * @param i
     *  row index to clear
     */
    virtual void clear(uint_t i);

    /**
     * Not made public, because the is_free method should be faster
     * @param i
     *  row index to check
     * @return
     *  true if entire row is zero
     */
    bool is_clear(uint_t i) const;

public:
    TableBase(uint_t row_size);

    TableBase(const TableBase& other);

    /**
     * the total number of records that can be stored in the BufferWindow.
     */
    uint_t capacity() const {
        return m_bw.m_nrow;
    }

    uint_t nrow_in_use() const {
        DEBUG_ASSERT_EQ(bool(m_bw.cend()), bool(m_bw.cbegin()),
                        "if buffer is set, the HWM should be set, else the HWM should not be set");
        return m_bw.cend() ? row_index(m_bw.cend()) : 0ul;
    }

    /**
     * clear the entire table only if it contains no protected records, else fatal error. this also resets the freed
     * rows objects, so it is made public
     */
    virtual void clear();

    /**
     * @return
     *  true if entire buffer window is zero (even above m_hwm)
     */
    bool is_clear() const;

    /**
     * strict equality: other table may have the same contents but a different buffer size, in this case the tables are
     * not considered equal. main use for this method is to check for successful copy operations
     * @param other
     *  table to compare against
     * @return
     *  true if the other table is identical in value to this one
     */
    bool operator==(const TableBase& other) const {
        if (this == &other) return true;
        return m_bw == other.m_bw;
    }

    bool operator!=(const TableBase& other) const {
        return !(*this==other);
    }
    /**
     * the size requirement of a single row in bytes (always an integer number of system words)
     */
    uint_t row_size() const {
        return m_bw.m_row_size;
    }

    bool i_can_modify() const {
        return m_bw.i_can_modify();
    }

    buf_t* trash_dst() {
        return m_bw.trash_dst();
    }

    /**
     * @return
     *  pointer to the first data word of the BufferWindow
     */
    buf_t *begin() {return m_bw.begin();}

    /**
     * @return
     * const pointer to the first data word of the BufferWindow
     */
    const buf_t *cbegin() const {return m_bw.cbegin();}

    /**
     * no end() method, since it is never deferencable
     * @return
     *  pointer to the first byte beyond the usable range of the BufferWindow
     */
    const buf_t *cend() const {return m_bw.cend();}

    /**
     * @param i
     *  row index
     * @return
     *  pointer to the beginning of the indexed row
     */
    buf_t *begin(uint_t i) {
        return m_bw.begin() + i * row_size();
    }

    /**
     * @param i
     *  row index
     * @return
     *  const pointer to the beginning of the indexed record
     */
    const buf_t *cbegin(uint_t i) const {
        return m_bw.cbegin() + i * row_size();
    }

    /**
     * Associate the table with a buffer by assigning the table an available BufferWindow
     * @param buffer
     *  pointer to buffer
     */
    void set_buffer(Buffer *buffer);

    /**
     * increase the high water mark, expand first if necessary
     * @param n
     *  number of records to bring into use
     * @return
     *  index of the first newly-accessible record
     */
    uint_t push_back(uint_t n=1ul);

    /**
     * If there are indices on the m_free_records stack: pop one and use it, else: push_back
     * @return
     *  index of free row
     */
    uint_t get_free_row();

    /**
     * clear the record of the indexed row and designate it freed. If it is already freed or protected, a fatal error is
     * thrown
     * @param i
     *  index of row to free
     */
    virtual void free(uint_t i);

    virtual void free();

    bool empty() const;

    const std::stack<uint_t>& freed_rows() const {
        return m_freed_rows;
    }

    uint_t nfreed_row() const {
        DEBUG_ASSERT_TRUE(freed_rows_consistent(), "inconsistent freed row indices");
        return m_freed_rows.size();
    }
    /**
     * the number of records (rows in use that have not been freed)
     * @return
     */
    uint_t nrecord() const {
        return nrow_in_use() - nfreed_row();
    }

    /**
     * @param i
     *  index of row
     * @return
     *  true if indexed row is freed
     */
    bool is_freed(uint_t i) const;

    /**
     * @return
     *  size of the buffer window in bytes
     */
    uint_t bw_size() const;

    /**
     * call the resize method on the buffer window and reflect the reallocation in m_nrow
     * @param nrow
     *  minimum number of rows in the new buffer.
     */
    virtual void resize(uint_t nrow, double factor=-1.0);

    /**
     * resize based on the number of additional rows required beyond those currently allocated
     * @param nrow
     *  minimum number of new rows
     */
    void expand(uint_t nrow, double factor=-1.0);

    /**
     * simply call free on each indexed row
     * @param irows
     *  all row indices marked for erasure
     */
    void free_many(const uintv_t& irows);

    /**
     * function pointer type for the callback associated with row transfers.
     * see RankAllocator.h
     */
    typedef std::function<void(const uintv_t&, uint_t, uint_t)> transfer_cb_t;
    /**
     * function pointer type for the callback associated with receipt of a single row in a transfer operation
     * see RankAllocator.h
     */
    typedef std::function<void(uint_t)> recv_cb_t;

#if 0
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
#endif

    /**
     * copy a single record from another "source" table
     * @param src
     *  source table which must have the same record size. should also have the same overall data layout but this is not
     *  verified at this level
     * @param isrc
     *  record in the source table to be copied
     * @param idst
     *  record in this table to which the source record is copied
     */
    void copy_record_in(const TableBase& src, uint_t isrc, uint_t idst);

    /**
     * swap record contents via std::swap_ranges
     * @param i
     *  index of record to be swapped
     * @param j
     *  index of other record to be swapped
     */
    void swap_records(uint_t i, uint_t j);

    /**
     * increase protection level of the record. insert key-value pair first if non-existent
     * @param i
     *  index of the record to protect
     */
    void protect(uint_t i) const {
        DEBUG_ASSERT_LT(i, nrow_in_use(), "row index OOB");
        DEBUG_ASSERT_FALSE(m_is_freed_row[i], "cannot protect a freed record");
        auto it = m_protected_rows.find(i);
        if (it == m_protected_rows.end()) m_protected_rows.insert({i, 1ul});
        else ++it->second;
    }

    /**
     * @param i
     *  record index
     * @return
     *  number of protect(i) calls - number of unprotect(i) calls
     */
    uint_t protection_level(uint_t i) const {
        DEBUG_ASSERT_LT(i, nrow_in_use(), "record index OOB");
        auto it = m_protected_rows.find(i);
        if (it == m_protected_rows.end()) return 0ul;
        else return it->second;
    }

    /**
     * @param i
     *  index of record
     * @return
     *  true if the protection level is non-zero
     */
    bool is_protected(uint_t i) const {
        return protection_level(i);
    }

    /**
     * decrement protection level of the indexed record
     * @param i
     *  index of protected record
     */
    void unprotect(uint_t i) const {
        DEBUG_ASSERT_TRUE(protection_level(i), "can't unprotect a record which already has zero protection");
        auto it = m_protected_rows.find(i);
        if (it->second==1ul) m_protected_rows.erase(it);
        else --it->second;
    }

    /**
     * debugging method to ensure consistency of the two structures keeping track of the free status of records
     * @return
     *  true is m_free_records member is consistent with m_is_free_record
     */
    bool freed_rows_consistent() const;

    /**
     * "Location" class which describes the location of a record in a distributed table i.e. by a local record index and
     * an MPI rank index
     */
    struct Loc {
        /**
         * rank and record indices
         */
        const uint_t m_irank, m_irec;

        Loc(uint_t irank, uint_t irec);

        /**
         * call on all procs, with all irecs = ~0ul except for on the owning rank
         * @param irec
         *  record (row) index in the distributed table if owned by MPI rank, else ~0ul
         */
        Loc(uint_t irec);

        /**
         * @return
         *  true if the location is anywhere i.e. has a valid rank index
         */
        operator bool() const;

        operator uintv_t() const {
            if (is_mine()) return {m_irec};
            return {};
        }

        /**
         * @return
         *  true if MPI rank index matches bcast-shared rank of identified record
         */
        bool is_mine() const;

        bool operator==(const Loc& other);

        bool operator!=(const Loc& other);
    };

    /**
     * When Field-based data structure is introduced in the derived classes, this method is capable of displaying
     * human-readable columns. here though, we don't know the data layout, so return record bytes as integers
     * @param ordering
     *  optional indices to reorder table on the fly as it is printed, without the need to physically reorder records
     * @return
     *  string representing table's contents
     */
    virtual str_t to_string(const uintv_t* ordering = nullptr) const;

    /**
     * gather contents of another table over all MPI ranks into this table on all ranks
     * @param src
     *  table whose records are to be gathered.
     */
    virtual void all_gatherv(const TableBase& src);

    /**
     * gather contents of another table over all MPI ranks into this table on the given rank
     * @param src
     *  table whose records are to be gathered.
     * @param irank
     *  index of the only rank in the communicator to receive the data from src
     */
    virtual void gatherv(const TableBase& src, uint_t irank = 0ul);

    /**
     * @return
     *  true if this Table contains any protected records
     */
    bool is_protected() const;

    str_t name() const {
        return m_bw.name();
    }

    void rename(str_t name) const {
        m_bw.rename(name);
    }

};


#endif //M7_TABLEBASE_H
