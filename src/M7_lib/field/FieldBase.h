//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_FIELDBASE_H
#define M7_FIELDBASE_H

#include <M7_lib/defs.h>
#include <cstring>
#include <utility>
#include <M7_lib/util/Hash.h>
#include <M7_lib/nd/NdArrayList.h>
#include <M7_lib/hdf5/NdDistList.h>
#include <M7_lib/hdf5/Node.h>

#include "Row.h"


/**
 * Base class for the basic containers of data within Rows, which in turn reference locations within a Buffer via Table.
 *
 * Due to the bi-directional references between Fields and Rows, and the dual role of these classes as data LAYOUT
 * specifiers and pointers to the DATA itself, special methods for copying inthis class and all derived classes must be
 * carefully implemented to respect the following semantics:
 *
 *  - copy ctor: the LAYOUT is being copied since the containing Row is delegating this ctor in the process of being
 *          itself copy-constructed. Raise an error if the Row's m_child member is undefined, and also copy assign so
 *          that the initial value of the copied Field is that of the Field from which it was copied
 *  - copy assign: the DATA is being copied. This method is not delegated by the Row, since Row data copies are done by
 *          buffer copy without reference to the layout defined by the associated Field set. As the main method of data
 *          copy, it is performance CRITICAL, therefore the base assigment definition and that of all subclasses should
 *          be defined in the header files to give the compiler the chance to inline
 *
 * move semantics are not utilized since all necessary use cases for Fields are covered by the rule-of-three methods.
 * The Field never needs to be moved into another symbol, since user-defined Rows and CompositeFields are classes with
 * their own members, and so the initializing ctor is sufficient.
 */
struct FieldBase {
    /**
     * Row object to which this field belongs, and has been apportioned a range of bytes within
     */
    Row *m_row = nullptr;
    /**
     * the derived classes which implement the features of particular types of Fields will carry type information, this
     * is passed in std::type_info form to the base class for logging purposes
     */
    const std::type_info &m_type_info;
    /**
     * number of bytes required to store the field
     * not necessarily an integral multiple of c_nbyte_word, although a field added to a row after another of a
     * different type will always begin on a word boundary
     */
    const uint_t m_size;
    /**
     * string identifier for logging purposes
     */
    const str_t m_name;
    /**
     * even if the next field added to a Row after this one is of the same type, it always begins at the system word
     * boundary in the buffer. this is to aid bytewise operations on this field without the need to be concerned with
     * encountering data corresponding to the next field in the Row
     */
    const bool m_force_own_words;
private:
    /**
     * a zero-valued string so that memory can be cleared by copy
     */
    v_t<buf_t> m_null_string;
    /**
     * number of bytes by which Field is offset from the beginning of the Row (first Field to be added has offset 0)
     */
    uint_t m_row_offset = ~0ul;
    /**
     * index within m_fields vector of m_row
     */
    uint_t m_row_index = ~0ul;
    friend Row;

public:

    /**
     * @param row
     *  row object to which the new field should be added
     * @param size
     *  bytes required to store the field (not including the padding between dissimilarly-typed fields in a Row)
     * @param type_info
     *  run time evaluated type signature, only used so the Row can determine whether to jump the offset to the next word
     * @param name
     *  name of the field i.e. "column heading"
     */
    FieldBase(Row *row, uint_t size, const std::type_info &type_info, str_t name, bool force_own_words);

    FieldBase(const FieldBase &other);

    FieldBase &operator=(const FieldBase &other){
        if (&other == this) return *this;
        DEBUG_ASSERT_TRUE(is_comparable(other), "can't compare to incompatible field");
        std::memcpy(begin(), other.begin(), m_size);
        return *this;
    }

    bool is_comparable(const FieldBase &other) const;

    bool belongs_to_row() const;

    bool belongs_to_row(const Row* row) const;

    bool belongs_to_row(const Row& row) const;

    buf_t *begin() const;

    buf_t *end() const;

    const Row *row() const;

    Row *row_of_copy() const;

    void zero();

    bool is_zero() const;

    bool operator==(const FieldBase &other) const;

    bool operator!=(const FieldBase &other) const;

    bool operator<(const FieldBase &other) const;

    bool operator>(const FieldBase &other) const;

    bool operator<=(const FieldBase &other) const;

    bool operator>=(const FieldBase &other) const;

    hash::digest_t hash() const;

    virtual str_t to_string() const = 0;

    uint_t to_buffer(buf_t* buf, uint_t irow_begin, uint_t nitem_max, std::set<uint_t> irows_empty) const;

    void to_buffer(buf_t* dst, uint_t irow) const {
        DEBUG_ASSERT_LT(irow, m_row->m_table->nrow_in_use(), "row index OOB");
        DEBUG_ASSERT_FALSE(m_row->m_table->is_freed(irow), "copying a freed row to buffer");
        auto src = m_row->m_table->begin() + m_row->m_size * irow + m_row_offset;
        std::memcpy(dst, src, m_size);
    }

    /**
     * save a entire field as raw bytes
     * @param nw
     *  NodeWriter instance (HDF5 group or file)
     * @param name
     *  desired name in the
     * @param this_rank
     *  true if data from this MPI rank is to be included in the saved set
     */
    virtual void save(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank) const;

    void save(const hdf5::NodeWriter& nw, bool this_rank) const;

    virtual uintv_t h5_shape() const {
        return {};
    }

    virtual strv_t h5_dim_names() const {
        return {};
    }

    virtual hdf5::Type h5_type() const {
        return {};
    }

private:

    int cmp(const FieldBase &other) const;
};

static std::ostream &operator<<(std::ostream &os, const FieldBase &v) {
    os << v.to_string();
    return os;
}

#endif //M7_FIELDBASE_H
