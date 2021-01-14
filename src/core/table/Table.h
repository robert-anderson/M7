//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <src/core/util/utils.h>
#include <src/core/field/Fields.h>
#include "src/defs.h"
#include "src/core/field/TableField.h"
#include "Buffer.h"
#include "src/core/field/Flag.h"
#include "src/core/sort/ExtremalValues.h"


struct Table {
    std::vector<const TableField *> m_fields;
    Buffer::Window m_bw;
    size_t m_row_size;
    size_t m_row_dsize;
    size_t m_current_byte_offset = 0ul;
    char *m_data;
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;

    void set_buffer(Buffer* buffer){
        ASSERT(buffer);
        ASSERT(!m_bw.dbegin())
        buffer->append_window(&m_bw);
    }

    bool is_full() const;

    size_t push_back(size_t nrow=1);



    defs::data_t* dbegin();

    const defs::data_t* dbegin() const;

    char *begin();

    const char *begin() const;



    defs::data_t* dbegin(const size_t &irow);

    const defs::data_t* dbegin(const size_t &irow) const;

    char *begin(const size_t &irow);

    const char *begin(const size_t &irow) const;

    size_t add_field(const TableField *field);

    void clear();

    void clear(const size_t &irow);

    size_t bw_dsize() const;

    std::string field_details(size_t width=30) const;

    void print_field_details(size_t width=30) const;

    void print_contents(const defs::inds* ordering=nullptr) const;

    void print_contents(const ExtremalValues& xv) const;

    bool is_cleared() const;

    bool is_cleared(const size_t &irow) const;

    virtual void erase_rows(const defs::inds& irows) {
        for (auto irow : irows) clear(irow);
    }

    void resize(size_t nrow){
        assert(nrow>m_nrow);
        m_bw.resize(nrow*m_row_dsize);
        m_nrow = nrow;
    }

    void expand(){
        m_bw.expand();
        m_nrow = m_bw.dsize()/m_row_dsize;
    }

    void expand(size_t nrow){
        assert(nrow>m_nrow);
        m_bw.expand(nrow*m_row_dsize);
        m_nrow +=nrow;
    }

#if 0
    virtual void insert_rows(const Buffer& recv){
        const auto nrow = recv.dsize()/m_row_dsize;
        for (size_t irow = 0; irow<nrow; ++irow){

        }
    }

    void transfer_rows(const defs::inds& irows, size_t irank_src, size_t irank_dst){
        size_t nrow;
        if (mpi::i_am(irank_src)) {
            nrow = irows.size();
            // make all rows to be sent contiguous in memory
            Buffer tmp("Outward transfer buffer", m_row_dsize, nrow);
            for (auto iirow = 0ul; iirow < nrow; ++iirow) {
                const auto &irow = irows[iirow];
                std::memcpy(tmp.ptr() + iirow * m_row_dsize, dbegin(irow), m_row_size);
            }
            mpi::send(&nrow, 1, irank_dst,0);
        }
        if (mpi::i_am(irank_dst)) {
            mpi::recv(&nrow, 1, irank_dst, 0);
            Buffer tmp("Inward transfer buffer", m_row_dsize, nrow);
            mpi::recv(tmp.ptr(), m_row_dsize * nrow, irank_dst, 1);
            // now emplace received rows in table buffer window
            for (auto iirow = 0ul; iirow < nrow; ++iirow) {
                const auto &irow = irows[iirow];
                std::memcpy(dbegin(irow), tmp.ptr() + iirow * m_row_dsize, m_row_size);
            }
        }
        // sent rows can now be removed
        if (mpi::i_am(irank_src)) {
            erase_rows(irows);
        }
    }
#endif
};


#endif //M7_TABLE_H
