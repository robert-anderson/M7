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
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;
    /*
     * when copying the Table, the fields being copied need to know
     * the correct m_table pointer, so we retain a pointer to the last
     * copied Table
     */
    mutable Table *m_last_copied = nullptr;

    Table();

    Table(const Table &other);

    void set_buffer(Buffer *buffer);

    bool is_full() const;

    size_t push_back(size_t nrow = 1);

    defs::data_t *dbegin();

    const defs::data_t *dbegin() const;

    char *begin();

    const char *begin() const;

    defs::data_t *dbegin(const size_t &irow);

    const defs::data_t *dbegin(const size_t &irow) const;

    char *begin(const size_t &irow);

    const char *begin(const size_t &irow) const;

    size_t add_field(const TableField *field);

    void clear();

    void clear(const size_t &irow);

    size_t bw_dsize() const;

    std::string field_details(size_t width = 30) const;

    void print_field_details(size_t width = 30) const;

    void print_contents(const defs::inds *ordering = nullptr) const;

    void print_contents(const ExtremalValues &xv) const;

    bool is_cleared() const;

    bool is_cleared(const size_t &irow) const;

    void resize(size_t nrow);

    void expand(size_t nrow);

    virtual void erase_rows(const defs::inds &irows);

    virtual void insert_rows(const Buffer &recv);

    void transfer_rows(const defs::inds &irows, size_t irank_src, size_t irank_dst);

    bool has_compatible_format(const Table &other);

    void copy_row(const Table &src, size_t irow_src, size_t irow_dst);

};


#endif //M7_TABLE_H
