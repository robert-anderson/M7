//
// Created by Robert J. Anderson on 08/05/2021.
//

#ifndef M7_TABLEBASE_H
#define M7_TABLEBASE_H

#include <M7_lib/io/Logging.h>

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
    const int m_nrow_p2p_tag = mpi::new_p2p_tag();
    const int m_irows_p2p_tag = mpi::new_p2p_tag();

    RowTransfer(str_t name) :
            m_send_buffer("Outward transfer buffer", 1),
            m_recv_buffer("Inward transfer buffer", 1) {
        logging::info("Initializing row send/recv buffers for table \"{}\"", name);
        logging::debug("P2P tag for number of row indices to transfer for \"{}\": {}", name, m_nrow_p2p_tag);
        logging::debug("P2P tag for array of row indices to transfer for \"{}\": {}", name, m_irows_p2p_tag);
        m_send_buffer.append_window(&m_send_bw);
        m_recv_buffer.append_window(&m_recv_bw);
    }
};

struct RowProtector;

/**
 * Base class for all row-contiguous data in the program.
 */
struct TableBase {
    /**
     * The data layout of a Table is defined by a class derived from Row, but this templated dependence is not included
     * in this base class. While the fine detail of field offsets etc is not needed here, we need to know the basics,
     * i.e. what amount of memory does a single row require. This is always padded to a whole number of system words.
     * A system word is sizeof(uint_t) bytes in length
     */
    /**
     * The contents of the Table are ultimately stored in a Buffer, which is just a wrapper for a dynamically
     * allocated array of data_t elements. Buffers can be shared by many Tables though, so a table is instead
     * given access to a "window" of that Buffer. All tables which have windows on the same buffer have the ability
     * to resize, and in response to this, the Buffer will ensure that table row data is moved to new positions
     * in the resized buffer, and that each of the other windows are pointed to the beginning of that data.
     *
     * the m_row_size member of m_bw is the definitive length of the Table's row in bytes, and the Row class ensures
     * that this is always an integer multiple of the system word length
     */
    Buffer::Window m_bw;
    /**
     * "high water mark" is result of the next call to push_back. Think of this as the Table's size by analogy to
     * v_t
     */
    uint_t m_hwm = 0ul;
    /**
     * indices of vacated rows below the high water mark should be pushed into this stack to allow reuse
     */
    std::stack<uint_t> m_free_rows;
    /**
     * buffered space for communication in the event of rank reallocation, instantiate on first transfer if required
     */
    std::unique_ptr<RowTransfer> m_transfer = nullptr;
    /**
     * save a buffer of zeros for fast memcmp to determine whether a row is cleared
     */
    const v_t<char> m_null_row_string;
    /**
     * if a row is to be protected, its index must appear as a key in the following map.
     * the value associated with the key is the "protection level" of the row. if one object requires that a row i not
     * be deleted, then {i: 1} appears in the map. If a second object requires that row i not be deleted, then the
     * record is looked up, then the value is incremented so than {i: 2} appears in the map. If then the first object
     * releases its protection of row i, the record returns to {i: 1}, so that the row will still not be deleted even
     * though the first dependent object has no objection. the row will remain protected until the second object
     * releases, and the record is deleted from the map
     */
    mutable std::map<uint_t, uint_t> m_protected_rows;


    TableBase(uint_t row_size);

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
        if (other.nrow()) {
            resize(other.nrow(), 0.0);
            m_bw = other.m_bw;
        }
        m_hwm = other.m_hwm;
        m_free_rows = other.m_free_rows;
        m_protected_rows = other.m_protected_rows;
        return *this;
    }

public:
    TableBase(const TableBase& other);

    /**
     * the number of rows in the BufferWindow. Think of this as the Table's capacity by analogy to v_t
     */
    uint_t nrow() const {
        return m_bw.m_nrow;
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
        if (nrow()!=other.nrow()) return false;
        if (m_hwm!=other.m_hwm) return false;
        return std::memcmp(begin(), other.begin(), m_hwm * row_size())==0;
    }

    bool operator!=(const TableBase& other) const {
        return !(*this==other);
    }
    /**
     * the size of a single row in bytes (always an integer number of system words)
     */
    uint_t row_size() const {
        return m_bw.m_row_size;
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
     * @param irow
     *  row index
     * @return
     *  pointer to the first data word of the indexed row in the BufferWindow
     */
    buf_t *begin(uint_t irow);

    /**
     * @param irow
     *  row index
     * @return
     *  const pointer to the first data word of the indexed row in the BufferWindow
     */
    const buf_t *begin(uint_t irow) const;

    /**
     * Associate the table with a buffer by assigning the table an available BufferWindow
     * @param buffer
     *  pointer to buffer
     */
    void set_buffer(Buffer *buffer);

    /**
     * @return
     *  true if the high water mark can't be increased without first increasing the number of rows
     */
    bool is_full() const;

    /**
     * increase the high water mark, expand first if necessary
     * @param nrow
     *  number of rows to make accessible
     * @return
     *  row index of the first newly-accessible row
     */
    uint_t push_back(uint_t nrow = 1);

    /**
     * If there are free rows on the stack: pop one and use it, else: push_back
     * @return
     *  index of unused row
     */
    uint_t get_free_row();

    /**
     * clear the entire table only if it contains no protected rows, else fatal error
     */
    virtual void clear();

    /**
     * clear the indexed row only if it is not protected, else fatal error
     * @param irow
     *  row index to clear
     */
    virtual void clear(uint_t irow);

    /**
     * @return
     *  size of the buffer window in bytes
     */
    uint_t bw_size() const;

    /**
     * @return
     *  true if entire buffer window is zero (even above m_hwm)
     */
    bool is_cleared() const;

    /**
     * @param irow
     *  row index to check
     * @return
     *  true if entire row is zero
     */
    bool is_cleared(uint_t irow) const;

    /**
     * call the resize method on the buffer window and reflect the reallocation in m_nrow
     * @param nrow
     *  minimum number of rows in the new buffer. the buffer's resize_factor determines the actual size of the reallocation
     */
    void resize(uint_t nrow, double factor=-1.0);

    /**
     * resize based on the number of additional rows required beyond those currently allocated
     * @param nrow
     *  minimum number of new rows. the buffer's resize_factor determines the actual size of the reallocation
     */
    void expand(uint_t nrow, double factor=-1.0);

    /**
     * simply call clear on each indexed row
     * @param irows
     *  all row indices marked for erasure
     */
    void clear_rows(const uintv_t& irows);

    /**
     * in some derived classes, there is more to adding new rows than simply copying their contents into the hwm. In
     * these cases we need to call a method after the copy in order to maintain the integrity of those data structures
     * see MappedTable.h for more details
     * @param iinsert
     *  row index of the already-inserted data
     */
    virtual void post_insert(uint_t /*iinsert*/){}

    /**
     * If we have copied a block of rows contiguously into a table with non-trivial post-insert obligations, we need to
     * call post_insert for each copied row
     * @param ibegin
     *  row index of beginning of the already-inserted data
     * @param iend
     *  row index of end of the already-inserted data
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
     * insert all rows held in the buffer window into this Table, then call all callbacks if any.
     * @param recv
     *  recieved buffer containing transferred rows contiguously, and densely - no zeroed rows.
     * @param nrow
     *  number of rows being inserted
     * @param callbacks
     *  functions to call each time a received row is processed
     */
    virtual void insert_rows(const Buffer::Window& recv, uint_t nrow, const std::list<recv_cb_t>& callbacks);

    /**
     * By P2P MPI communication, send the rows identified in the first arg from irank_send to irank_recv.
     * Called on all ranks, but the first arg is only respected on the irank_send rank.
     * @param irows
     *  row indices on irank_send to be found, copied to contiguous buffer and transferred
     * @param irank_send
     *  MPI rank index from which the rows are being transferred
     * @param irank_recv
     *  MPI rank index to which the rows are being transferred
     * @param callbacks
     *  functions to call on irank_recv each time a received row is processed
     */
    void transfer_rows(const uintv_t& irows, uint_t irank_send, uint_t irank_recv,
                       const std::list<recv_cb_t>& callbacks = {});

    /**
     * copy a single row from another "source" table
     * @param src
     *  source table which must have the same row length. should also have the same overall data layout but this is not
     *  verified at this level
     * @param irow_src
     *  row in the source table to be copied
     * @param irow_dst
     *  row in this table to which the source row is copied bytewise
     */
    void copy_row_in(const TableBase& src, uint_t irow_src, uint_t irow_dst);

    /**
     * swap row contents dword-for-dword via std::swap
     * @param irow
     *  row index to be swapped
     * @param jrow
     *  row index to be swapped
     */
    void swap_rows(uint_t irow, uint_t jrow);


    void protect(uint_t irow) const {
        auto it = m_protected_rows.find(irow);
        if (it==m_protected_rows.end()) m_protected_rows.insert({irow, 1ul});
        else ++it->second;
    }

    uint_t protection_level(uint_t irow) const {
        auto it = m_protected_rows.find(irow);
        if (it==m_protected_rows.end()) return 0ul;
        else return it->second;
    }

    bool is_protected(uint_t irow) const {
        return protection_level(irow);
    }

    void release(uint_t irow) const {
        DEBUG_ASSERT_TRUE(protection_level(irow), "can't release an unprotected row");
        auto it = m_protected_rows.find(irow);
        if (it->second==1ul) m_protected_rows.erase(it);
        --it->second;
    }


    /**
     * "Location" class which describes the location of a row in a distributed table i.e. by a row index and a rank index
     */
    struct Loc {
        /**
         * rank and row indices
         */
        const uint_t m_irank, m_irow;

        Loc(uint_t irank, uint_t irow);

        /**
         * @return
         *  true if the location is anywhere i.e. has a valid rank index
         */
        operator bool() const;

        /**
         * @return
         *  true if MPI rank index matches bcast-shared rank of identified row
         */
        bool is_mine() const;

        bool operator==(const Loc& other);

        bool operator!=(const Loc& other);
    };

    /**
     * When Field-based data structure is introduced in the derived classes, this method is capable of displaying
     * human-readable columns. here though, we don't know the data layout, so return bytes as integers
     * @param ordering
     *  optional indices to reorder table on the fly as it is printed, without the need to physically reorder rows
     * @return
     *  string representing table's contents
     */
    virtual str_t to_string(const uintv_t* ordering = nullptr) const;

    /**
     * gather contents of another table over all MPI ranks into this table on all ranks
     * @param src
     *  table whose rows are to be gathered.
     */
    virtual void all_gatherv(const TableBase& src);

    /**
     * gather contents of another table over all MPI ranks into this table on the given rank
     * @param src
     *  table whose rows are to be gathered.
     * @param irank
     *  index of the only rank in the communicator to receive the data from src
     */
    virtual void gatherv(const TableBase& src, uint_t irank = 0ul);

    /**
     * @return
     *  true if any of the associated RowProtectors are protecting any rows of this table
     */
    bool is_protected() const;

    /**
     * @return
     *  number rows below the high water mark that aren't free
     */
    uint_t nrow_nonzero() const;

    str_t name() const {
        return m_bw.name();
    }
};


#endif //M7_TABLEBASE_H
