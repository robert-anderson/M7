//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_ROW_H
#define M7_ROW_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/table/TableBase.h>
#include <M7_lib/util/Pointer.h>

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
     * the row needs access to the buffer window and high water mark of the table with which it is associated
     */
    TableBase *m_table;
    /**
     * vector of the associated fields. These are polymorphic pointers to the symbols defined in the subclasses
     */
    v_t<FieldBase *> m_fields;
private:
    /**
     * the position of the currently selected record in the
     */
    mutable buf_t *m_begin = nullptr;

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
     * the m_begin pointer must be in the range [table begin, table hwm) for valid dereferencing (field access)
     * @return
     *  true if the m_begin pointer is valid with respect to the "in use" range of the table object
     */
    bool is_deref_valid() const {
        if (!m_table->m_bw.m_begin) return false;
        DEBUG_ASSERT_TRUE(m_table->m_hwm, "buffer window has null HWM pointer");
        return ptr::in_range(m_begin, m_table->m_bw.m_begin, m_table->m_hwm);
    }

    /**
     * the Row can be in a state valid for record indexing, even when invalid for dereferencing. this is a debugging
     * method used to ensure the Row is in a valid state
     * @return
     *  true if the m_begin pointer is dereferencable or (it is nullptr and so is table begin) or (it is equal to the
     *  high water mark) the latter condition is typically fulfilled at the end of a loop over rows
     */
    bool is_valid() const {
        return (m_begin==m_table->m_bw.m_begin) || is_deref_valid() || m_begin == m_table->m_hwm;
    }

    operator bool () const {
        return is_deref_valid();
    }

    /**
     * @return
     * record position within Table
     */
    uint_t index() const {
        DEBUG_ASSERT_TRUE(is_valid(), "the row is not pointing to memory in the permitted range");
        const auto begin_offset = std::distance(m_table->m_bw.m_begin, m_begin);
        return begin_offset / m_size;
    }


    /**
     * @param irow_end
     *  exclusive maximum value for the stored row index
     * @return
     *  true if stored row index is within the range [0, irow_end)
     */
    bool in_range(uint_t irow_end) const {
        DEBUG_ASSERT_TRUE(is_valid(), "invalid range end given");
        return ptr::in_range(m_begin, m_table->m_bw.m_begin, m_table->begin(irow_end));
    }

    buf_t *begin() {
        DEBUG_ASSERT_TRUE(m_begin, "the row pointer is not set")
        DEBUG_ASSERT_TRUE(is_deref_valid(), "the row is not pointing to memory in the dereferencable range");
        return m_begin;
    }

    const buf_t *begin() const {
        DEBUG_ASSERT_TRUE(m_begin, "the row pointer is not set")
        DEBUG_ASSERT_TRUE(is_deref_valid(), "the row is not pointing to memory in the dereferencable range");
        return m_begin;
    }

    /*
     * the 3 "cursor" methods
     */
    void restart(uint_t irow_begin) const {
        DEBUG_ASSERT_LE(irow_begin, m_table->nrow_in_use(), "Cannot restart to an out-of-range row index");
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        if (!m_table->m_hwm && !irow_begin){
            m_begin = nullptr;
        } else {
            DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
            m_begin = m_table->begin(irow_begin);
        }
    }

    void restart() const {
        restart(0);
    }

    /**
     * prefix increment to advance row to the next record
     */
    const Row& operator ++() const {
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
        DEBUG_ASSERT_TRUE(is_valid(), "Row is out of table bounds");
        m_begin += m_size;
        return *this;
    }

    void jump(uint_t i) const {
        DEBUG_ASSERT_TRUE(m_table, "Row must be assigned to a Table");
        DEBUG_ASSERT_TRUE(m_table->begin(), "Row is assigned to Table buffer window without a beginning");
        DEBUG_ASSERT_LE(i, m_table->nrow_in_use(), "Row is out of table bounds");
        m_begin = m_table->begin(i);
    }

    void jump(const Row &other) const {
        jump(other.index());
    }

    void push_back_jump() {
        jump(m_table->push_back());
    }

    void select_null() const {
        m_begin = nullptr;
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

    /**
     * increment protection level of the currently indexed row
     */
    void protect();
    /**
     * @return
     *  (number of times protect() has been called) - (number of times release() has been called) on the currently
     *  indexed row
     */
    uint_t protection_level() const;
    bool is_protected() const;
    /**
     * decrement protection level of the currently indexed row
     */
    void unprotect();

    str_t field_names_string() const;

    str_t to_string() const;

    uint_t add_field(FieldBase *field);

    uint_t nfield() const;

    void free();

    bool is_freed() const {
        return m_table->is_freed(index());
    }

    virtual bool is_h5_write_exempt() const;

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
    SingleFieldRow(Args... args): m_field(this, args...){}

    SingleFieldRow(const SingleFieldRow& other): Row(other), m_field(other.m_field){}

    field_t &key_field() {
        return m_field;
    };
};


#endif //M7_ROW_H
