//
// Created by rja on 09/02/2021.
//

#ifndef M7_ROW_H
#define M7_ROW_H

#include "src/core/util/utils.h"
#include "src/core/io/HDF5Wrapper.h"
#include "src/core/parallel/MPIAssert.h"
#include "src/core/table/TableBase.h"

struct FieldBase;

/**
 * Classes derived from Row define the data layout of Tables. Row-derived classes declare members derived from FieldBase
 * which give numerical meaning to the raw data stored in the buffer window assigned to a Table. Rows have a dual
 * purpose:
 *  1. assign offsets to the FieldBase-derived members
 *  2. act as a cursor pointing to a specific row in a buffer window
 *
 * The latter usage requires care. There are some potential ways this usage could have been avoided, but they are
 * problematic.
 *  1. have a table be represented by a vector of row objects, with each storing their own index: this would incur the
 *     cost in memory of at least an additional reference per table row, also the cost of extending table buffers would
 *     increase noticeably
 *  2. have a RowDefinition class and a RowCursor class, where RowDefinition does not store positional information, and
 *     RowCursor only requires a reference to the definition and a row index, so it can be constructed cheaply. The
 *     problem is that there is no way for the field objects stored in the RowDefinition to access this positional
 *     information.
 *
 * Overall, the statefulness of Rows is not a problem provided the programmer is aware which rows are in use in which
 * contexts. More often than not, the row object will be used at the beginning of a long loop, in which case the cost of
 * copying the row is negligible in comparison.
 */
struct Row {
    /**
     * the row needs access to the buffer window and high water mark of the table with which it is asscociated
     */
    TableBase *m_table;
    /**
     * these fields define the position of the currently selected row. when set to these initial "null" values, the
     * row is not pointing to anything and field dereferences / calls to step should fail asserts in the debug builds.
     * no checks are enforced in the release build.
     * TODO: investigate whether dbegin needs to be cached with regard to performance
     */
    mutable defs::data_t *m_dbegin = nullptr;
    mutable size_t m_i = ~0ul;
    /**
     * vector of the associated fields. These are polymorphic pointers to the symbols defined in the subclasses
     */
    std::vector<FieldBase *> m_fields;
    /**
     * total size of a single row in bytes
     */
    size_t m_size = 0ul;
    /**
     * total size of a single row in data words (defs::data_t)
     */
    size_t m_dsize = 0ul;
    /**
     * number of bytes already allocated to added fields
     */
    size_t m_current_offset = 0ul;
    /**
     * required only when copying
     */
    mutable Row *m_child = nullptr;

    /**
     * @param irow_end
     *  exclusive maximum value for the stored row index
     * @return
     *  true if stored row index is within the range [0, irow_end)
     */
    bool in_range(const size_t &irow_end) const {
        return m_i < irow_end;
    }

    /**
     * @return
     *  true if the stored row index is within the range [0, m_table->m_hwm)
     */
    bool in_range() const {
        return m_i < m_table->m_hwm;
    }

    /**
     * @return
     *  true if the m_dbegin pointer is valid with respect to the table object
     */
    bool ptr_in_range() const {
        return (m_dbegin >= m_table->dbegin()) && (m_dbegin < m_table->dbegin() + m_table->m_hwm * m_dsize);
    }

    defs::data_t *dbegin() {
        MPI_ASSERT(m_dbegin, "the row pointer is not set")
        MPI_ASSERT(ptr_in_range(), "the row is not pointing to memory in the permitted range");
        return m_dbegin;
    }

    const defs::data_t *dbegin() const {
        MPI_ASSERT(m_dbegin, "the row pointer is not set")
        MPI_ASSERT(ptr_in_range(), "the row is not pointing to memory in the permitted range");
        return m_dbegin;
    }

    /*
     * the 3 "cursor" methods
     */
    void restart(const size_t &irow_begin) const {
        MPI_ASSERT(irow_begin < m_table->m_hwm, "Cannot restart to an out-of-range row index");
        MPI_ASSERT(m_table, "Row must be assigned to a Table");
        MPI_ASSERT(m_table->dbegin(), "Row is assigned to Table buffer window without a beginning");
        m_dbegin = m_table->dbegin(irow_begin);
        m_i = irow_begin;
    }

    void restart() const {
        restart(0);
    }

    void step() const {
        MPI_ASSERT(m_table, "Row must be assigned to a Table");
        MPI_ASSERT(m_table->dbegin(), "Row is assigned to Table buffer window without a beginning");
        MPI_ASSERT(in_range(), "Row is out of table bounds");
        m_dbegin += m_dsize;
        m_i++;
    }

    void jump(const size_t &i) const {
        MPI_ASSERT(m_table, "Row must be assigned to a Table");
        MPI_ASSERT(m_table->dbegin(), "Row is assigned to Table buffer window without a beginning");
        m_dbegin = m_table->dbegin() + m_dsize * i;
        m_i = i;
        MPI_ASSERT(in_range(), "Row is out of table bounds");
    }

    void jump(const Row &other) const {
        jump(other.m_i);
    }

    void push_back_jump() {
        jump(m_table->push_back());
    }

    /**
     * if m_table has been reallocated, the cached m_dbegin pointer is not valid, calling this method rectifies that
     */
    void refresh() const {
        jump(m_i);
    }

    void copy_in(const Row &other) {
        ASSERT(other.m_dsize == m_dsize);
        std::copy(other.dbegin(), other.dbegin() + m_dsize, dbegin());
    }

    Row() {}

    Row(const Row &other) {
        m_table = other.m_table;
        m_dbegin = other.m_dbegin;
        other.m_child = this;
        ASSERT(m_fields.empty())
    }

    Row &operator=(const Row &other) {
        copy_in(other);
        return *this;
    }

    std::string field_names_string() const;

    std::string to_string() const;

    size_t add_field(FieldBase *field);

    void clear();

    bool is_cleared() const;

    virtual bool is_h5_write_exempt() const;

    bool is_protected() const;

};

template<typename row_t>
struct KeyField {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    typedef typename std::remove_reference<typename std::result_of<decltype(&row_t::key_field)(
            row_t)>::type>::type type;

    static type &get(row_t &row) { return row.key_field(); }

    static const type &get(const row_t &row) { return const_cast<row_t &>(row).key_field(); }
};


#endif //M7_ROW_H
