//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <src/core/util/utils.h>
#include <stack>
#include "src/defs.h"
#include "Buffer.h"
#include "src/core/field/Column.h"
#include "src/core/field/Flag.h"
#include "src/core/sort/ExtremalValues.h"


struct RowTransfer {
    // make rows to be sent contiguous in memory
    Buffer m_send_buffer, m_recv_buffer;
    Buffer::Window m_send_bw, m_recv_bw;
    const int m_nrow_p2p_tag = mpi::new_p2p_tag();
    const int m_irows_p2p_tag = mpi::new_p2p_tag();
    RowTransfer(std::string name):
    m_send_buffer("Outward transfer buffer", 1),
    m_recv_buffer("Inward transfer buffer", 1) {
        log::info("Initializing row send/recv buffers for table \"{}\"", name);
        log::debug("P2P tag for number of row indices to transfer for \"{}\": {}", name, m_nrow_p2p_tag);
        log::debug("P2P tag for array of row indices to transfer for \"{}\": {}", name, m_irows_p2p_tag);
        m_send_buffer.append_window(&m_send_bw);
        m_recv_buffer.append_window(&m_recv_bw);
    }
};

struct Table {
    std::vector<const ColumnBase *> m_columns;
    Buffer::Window m_bw;
    size_t m_row_size;
    size_t m_row_dsize;
    size_t m_current_byte_offset = 0ul;
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;
    /*
     * indices of vacated rows below the high water mark should be pushed
     * into this stack to allow reuse
     */
    std::stack<size_t> m_free_rows;
    /*
     * when copying the Table, the columns being copied need to know
     * the correct m_table pointer, so we retain a pointer to the last
     * copied Table
     */
    mutable Table *m_last_copied = nullptr;
    // instantiate on first transfer if required
    std::unique_ptr<RowTransfer> m_transfer = nullptr;

    Table();

    Table(const Table &other);

    void set_buffer(Buffer *buffer);

    bool is_full() const;

    size_t push_back(size_t nrow = 1);

    size_t get_free_row();

    defs::data_t *dbegin();

    const defs::data_t *dbegin() const;

    defs::data_t *dbegin(const size_t &irow);

    const defs::data_t *dbegin(const size_t &irow) const;

    size_t add_column(const ColumnBase *column);

    void clear();

    void clear(const size_t &irow);

    size_t bw_dsize() const;

    std::string column_details(size_t width = 30) const;

    void print_column_details(size_t width = 30) const;

    void print_contents(const defs::inds *ordering = nullptr) const;

    void print_contents(const ExtremalValues &xv) const;

    bool is_cleared() const;

    bool is_cleared(const size_t &irow) const;

    void resize(size_t nrow);

    void expand(size_t nrow, double expansion_factor);

    void expand(size_t nrow);

    virtual void erase_rows(const defs::inds &irows);

    virtual void post_insert(const size_t& iinsert);

    typedef std::function<void(const defs::inds&, size_t, size_t)> transfer_cb_t;
    typedef std::function<void(size_t)> recv_cb_t;

    virtual void insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks);

    void transfer_rows(const defs::inds &irows, size_t irank_send, size_t irank_recv, const std::list<recv_cb_t> & callbacks= {});

    bool has_compatible_format(const Table &other);

    void copy_row_in(const Table &src, size_t irow_src, size_t irow_dst);

    struct Loc {
        const size_t m_irank, m_irow;
        Loc(size_t irank, size_t irow);
        operator bool() const;
        bool is_mine() const;
        bool operator==(const Loc& other);
        bool operator!=(const Loc& other);
    };

};


#endif //M7_TABLE_H
