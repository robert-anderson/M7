//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_ROW_H
#define M7_ROW_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/table/TableBase.h>

struct FieldBase;

/**
 * Classes derived from Row define the data layout of Tables. Row-derived classes declare members derived from FieldBase
 * which give numerical meaning to the raw data stored in the buffer window assigned to a Table. Rows have a dual
 * purpose:
 *  1. assign byte offsets to the FieldBase-derived members
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
 * Overall, the statefulness of Rows is not a problem provided the programmer is aware which Row instances are in use in
 * which contexts. More often than not, the row object will be used at the beginning of a long loop, in which case the
 * cost of copying the row is negligible in comparison.
 */
struct Row {
    /**
     * the row needs access to the buffer window and high water mark of the table with which it is asscociated
     */
    TableBase *m_table;
    /**
     * vector of the associated fields. These are polymorphic pointers to the symbols defined in the subclasses
     */
    v_t<FieldBase *> m_fields;
private:
    /**
     * these fields define the position of the currently selected row. when set to these initial "null" values, the
     * row is not pointing to anything and field dereferences / calls to step should cause ASSERTs to fail in the debug
     * builds. no checks are enforced in the release build.
     *
     * if m_i < m_table->m_hwm:
     *  m_begin should not be nullptr
     * else if m_i == m_table->m_hwm:
     *  m_begin should be nullptr, and field dereferencing should fail ASSERTs
     * else:
     *  this state is invalid
     *
     * TODO: investigate whether dbegin needs to be cached with regard to performance
     */
    mutable buf_t *m_begin = nullptr;
    mutable uint_t m_i = 0ul;

public:
    /**
     * total size of a single row in bytes
     */
    uint_t m_size = 0ul;
    /**
     * number of bytes already allocated to added fields
     */
    uint_t m_current_offset = 0ul;
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
    bool in_range(const uint_t &irow_end) const {
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
        return (m_begin >= m_table->begin()) && (m_begin < m_table->begin() + m_table->m_hwm * m_size);
    }

    buf_t *begin() {
        DEBUG_ASSERT_TRUE(m_begin, "the row pointer is not set")
        DEBUG_ASSERT_TRUE(ptr_in_range(), "the row is not pointing to memory in the permitted range");
        return m_begin;
    }

    const buf_t *begin() const {
        DEBUG_ASSERT_TRUE(m_begin, "the row pointer is not set")
        DEBUG_ASSERT_TRUE(ptr_in_range(), "the row is not pointing to memory in the permitted range");
        return m_begin;
    }

    /**
     * m_i == m_table->m_hwm is not valid for access, but is the state in which the row position data are left at the
     * loop when the in_range() loop termination condition becomes false, so the assert doesn't fail in this case
     * @return
     *  row position within Table
     */
    const uint_t& index() const {
        DEBUG_ASSERT_LE(m_i, m_table->m_hwm, "the row index is not in the permitted range");
        return m_i;
    }

    /*
     * the 3 "cursor" methods
     */
    void restart(const uint_t &irow_begin) const {
        DEBUG_ASSERT_LE(irow_begin, m_table->m_hwm, "Cannot restart to an out-of-range row index");
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        if (!m_table->m_hwm && !irow_begin){
            m_begin = nullptr;
        } else {
            DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
            m_begin = m_table->begin(irow_begin);
        }
        m_i = irow_begin;
    }

    void restart() const {
        restart(0);
    }

    void step() const {
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
        DEBUG_ASSERT_TRUE(in_range(), "Row is out of table bounds");
        m_begin += m_size;
        m_i++;
    }

    void jump(const uint_t &i) const {
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
        m_begin = m_table->begin() + m_size * i;
        m_i = i;
        DEBUG_ASSERT_LE(i, m_table->m_hwm, "Row is out of table bounds");
    }

    void jump(const Row &other) const {
        jump(other.m_i);
    }

    void push_back_jump() {
        jump(m_table->push_back());
    }

    void select_null() const {
        m_i = m_table->m_hwm;
        m_begin = nullptr;
    }

    bool null_selected() const {
        DEBUG_ASSERT_EQ(m_i==m_table->m_hwm, m_begin==nullptr, "Row in inconsistent state");
        return m_begin;
    }

    /**
     * if m_table has been reallocated, the cached m_dbegin pointer is not valid, calling this method rectifies that
     */
    void refresh() const {
        jump(m_i);
    }

    void copy_in(const Row &other) {
        ASSERT(other.m_size == m_size);
        std::copy(other.begin(), other.begin() + m_size, begin());
    }

    Row() {}

    ~Row();

    Row(const Row &other);

    Row &operator=(const Row &other) {
        copy_in(other);
        return *this;
    }

    str_t field_names_string() const;

    str_t to_string() const;

    uint_t add_field(FieldBase *field);

    uint_t nfield() const;

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


template<typename field_t>
struct SingleFieldRow : Row {
    field_t m_field;
    template<typename ...Args>
    SingleFieldRow(Args... args): Row(), m_field(this, args...){}

    field_t &key_field() {
        return m_field;
    };
};

#endif //M7_ROW_H
