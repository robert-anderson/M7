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

    typedef std::list<std::function<void(const size_t&)>> cb_list_t;

    virtual void erase_rows(const defs::inds &irows, const cb_list_t& callbacks);

    virtual void insert_rows(const Buffer &recv, const cb_list_t& callbacks);

    void send_rows(const defs::inds &irows, size_t irank_dst, const cb_list_t& callbacks={});

    void recv_rows(const size_t irank_src, const cb_list_t& callbacks={});

    bool has_compatible_format(const Table &other);

    void copy_row_in(const Table &src, size_t irow_src, size_t irow_dst);

    struct Loc {
        const size_t m_irank, m_irow;
        Loc(size_t irank, size_t irow): m_irank(irank), m_irow(irow){}

        operator bool() const {
            return m_irank!=~0ul;
        }
        bool is_mine() const {
            return mpi::i_am(m_irank);
        }
        bool operator==(const Loc& other){
            return m_irank==other.m_irank and m_irow==other.m_irow;
        }
        bool operator!=(const Loc& other){
            return !(m_irank==other.m_irank);
        }
    };

};


#endif //M7_TABLE_H
