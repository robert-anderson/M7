//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLEBASEX_H
#define M7_TABLEBASEX_H

#include <src/core/util/utils.h>
#include <stack>
#include <src/core/field/FieldX.h>
#include <src/core/sort/ExtremalValues.h>
#include "src/defs.h"
#include "Buffer.h"
#include "src/core/io/Logging.h"

#if 0
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
#endif

struct TableBaseX {
    const defs::inds m_field_format;
    const size_t m_row_dsize;
    const size_t m_row_size;
    Buffer::Window m_bw;
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
    // instantiate on first transfer if required
    std::unique_ptr<RowTransfer> m_transfer = nullptr;

    TableBaseX(defs::inds field_format, size_t row_dsize);

    TableBaseX(const TableBaseX &other);

    defs::data_t *dbegin() {
        return m_bw.m_dbegin;
    }

    const defs::data_t *dbegin() const {
        return m_bw.m_dbegin;
    }

    defs::data_t *dbegin(const size_t &irow) {
        return m_bw.m_dbegin + irow * m_row_dsize;
    }

    const defs::data_t *dbegin(const size_t &irow) const {
        return m_bw.m_dbegin + irow * m_row_dsize;
    }


    void set_buffer(Buffer *buffer);

    bool is_full() const;

    size_t push_back(size_t nrow = 1);

    size_t get_free_row();

    void clear();

    void clear(const size_t &irow);

    size_t bw_dsize() const;

    //void print_contents(const defs::inds *ordering = nullptr) const;

    //void print_contents(const ExtremalValues &xv) const;

    bool is_cleared() const;

    bool is_cleared(const size_t &irow) const;

    void resize(size_t nrow);

    void expand(size_t nrow, double expansion_factor);

    void expand(size_t nrow);

    virtual void erase_rows(const defs::inds &irows);

    virtual void post_insert(const size_t &iinsert);

    typedef std::function<void(const defs::inds &, size_t, size_t)> transfer_cb_t;
    typedef std::function<void(size_t)> recv_cb_t;

    virtual void insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks);

    void transfer_rows(const defs::inds &irows, size_t irank_send, size_t irank_recv,
                       const std::list<recv_cb_t> &callbacks = {});

    bool has_compatible_format(const TableBaseX &other);

    void copy_row_in(const TableBaseX &src, size_t irow_src, size_t irow_dst);

    struct Loc {
        const size_t m_irank, m_irow;

        Loc(size_t irank, size_t irow);

        operator bool() const;

        bool is_mine() const;

        bool operator==(const Loc &other);

        bool operator!=(const Loc &other);
    };
};

template<typename row_t>
struct TableX : TableBaseX {
    static_assert(std::is_base_of<RowX, row_t>::value, "Template arg must be derived from Row");
    row_t m_row;

    TableX(row_t row) :
            TableBaseX(static_cast<const RowX &>(row).field_format(),
                       static_cast<const RowX &>(row).m_dsize), m_row(row) {
        m_row.set_table(this);
    }


    void print_contents(const defs::inds *ordering= nullptr) const {
        const auto n = ordering ? std::min(ordering->size(), m_hwm) : m_hwm;
        m_row.select_first();
        for (size_t iirow = 0ul; iirow < n; ++iirow) {
            auto irow = ordering ? (*ordering)[iirow] : iirow;
            std::cout << irow << ". " << m_row.to_string() << "\n";
            m_row.select_next();
        }
        std::cout << std::endl;
    }

    void print_contents(const ExtremalValues &xv) const {
        defs::inds tmp;
        tmp.reserve(xv.nfound());
        for (size_t i = 0ul; i < xv.nfound(); ++i) tmp.push_back(xv[i]);
        print_contents(&tmp);
    }

private:
    RowX & base_row() {
        return static_cast<RowX &>(m_row);
    }
};


#endif //M7_TableX_H
