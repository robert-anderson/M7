//
// Created by rja on 08/05/2021.
//

#ifndef M7_TABLEBASE_H
#define M7_TABLEBASE_H

#include "src/core/table/Buffer.h"
#include "src/core/io/Logging.h"

struct RowTransfer {
    // make rows to be sent contiguous in memory
    Buffer m_send_buffer, m_recv_buffer;
    Buffer::Window m_send_bw, m_recv_bw;
    const int m_nrow_p2p_tag = mpi::new_p2p_tag();
    const int m_irows_p2p_tag = mpi::new_p2p_tag();

    RowTransfer(std::string name) :
            m_send_buffer("Outward transfer buffer", 1),
            m_recv_buffer("Inward transfer buffer", 1) {
        log::info("Initializing row send/recv buffers for table \"{}\"", name);
        log::debug("P2P tag for number of row indices to transfer for \"{}\": {}", name, m_nrow_p2p_tag);
        log::debug("P2P tag for array of row indices to transfer for \"{}\": {}", name, m_irows_p2p_tag);
        m_send_buffer.append_window(&m_send_bw);
        m_recv_buffer.append_window(&m_recv_bw);
    }
};

struct TableBase {
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

    TableBase(size_t row_dsize);

    TableBase(const TableBase &other);

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

    virtual void clear();

    virtual void clear(const size_t &irow);

    size_t bw_dsize() const;

    bool is_cleared() const;

    bool is_cleared(const size_t &irow) const;

    void resize(size_t nrow);

    void expand(size_t nrow, double expansion_factor);

    void expand(size_t nrow);

    virtual void erase_rows(const defs::inds &irows);

    virtual void post_insert(const size_t &iinsert);

    void post_insert_range(size_t ibegin=0ul, size_t iend=~0ul){
        if (iend==~0ul) iend = m_hwm;
        for (size_t i=ibegin; i<iend; ++i) post_insert(i);
    }

    typedef std::function<void(const defs::inds &, size_t, size_t)> transfer_cb_t;
    typedef std::function<void(size_t)> recv_cb_t;

    virtual void insert_rows(const Buffer::Window &recv, size_t nrow, const std::list<recv_cb_t> &callbacks);

    void transfer_rows(const defs::inds &irows, size_t irank_send, size_t irank_recv,
                       const std::list<recv_cb_t> &callbacks = {});

    void copy_row_in(const TableBase &src, size_t irow_src, size_t irow_dst);

    void swap_rows(const size_t &irow, const size_t &jrow);

    struct Loc {
        const size_t m_irank, m_irow;

        Loc(size_t irank, size_t irow);

        operator bool() const;

        bool is_mine() const;

        bool operator==(const Loc &other);

        bool operator!=(const Loc &other);
    };

    virtual std::string to_string(const defs::inds *ordering = nullptr) const {
        return "";
    }

    virtual void all_gatherv(const TableBase& src) {
        defs::inds nrows(mpi::nrank());
        defs::inds counts(mpi::nrank());
        defs::inds displs(mpi::nrank());
        ASSERT(src.m_row_dsize==m_row_dsize);
        mpi::all_gather(src.m_hwm, nrows);
        counts = nrows;
        for (auto& v: counts) v*=m_row_dsize;
        mpi::counts_to_displs_consec(counts, displs);
        auto nrow_total = std::accumulate(nrows.cbegin(), nrows.cend(), 0ul);
        resize(nrow_total);
        mpi::all_gatherv(src.dbegin(), m_hwm*m_row_dsize, dbegin(), counts, displs);
    }
};


#endif //M7_TABLEBASE_H
