//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDBASE_H
#define M7_FIELDBASE_H

#include <src/defs.h>
#include <cstring>
#include <src/core/hash/Hashing.h>
#include <src/core/nd/NdArrayList.h>
#include "src/core/io/HDF5Wrapper.h"
#include "Row.h"

/**
 * Base class for the basic containers of data within Rows, which in turn reference locations within a Buffer via Table.
 *
 * Due to the bi-directional references between Fields and Rows, and the dual role of these classes as data-LAYOUT
 * specifiers and pointers to the DATA itself, copy/move special methods for this class and all derived classes must be
 * carefully implemented to respect the following semantics
 *
 *  - copy ctor: the LAYOUT is being copied since the containing Row is delegating this ctor in the process of being
 *          itself copy-constructed. Raise an error if the Row If the Row's m_child member is undefined
 *  - copy assign: the DATA is being copied. This method is not delegated by the Row, since Row data copies are done by
 *          buffer copy without reference to the layout defined by the associated Field set
 *  - move ctor: the LAYOUT is being moved into a new symbol. Here, the row pointer is copied from the FieldBase being
 *          moved since it is assumed that this moving is due to forwarding in the execution of Row and CompoundField
 *          ctors.
 *  - move assign: semantically, the DATA is being moved, but the Field doesn't own the data buffer, it only points
 *          to a certain location within it. Attempted move assignment of Fields could be considered an error, but it is
 *          fine to simply delegate to the copy assignment instead.
 */
struct FieldBase {
    Row *m_row = nullptr;
    const std::type_info &m_type_info;
    const size_t m_size;
    const std::string m_name;
    size_t m_row_offset = ~0ul;
    std::vector<char> m_null_string;

private:
    size_t m_row_index = ~0ul;

    FieldBase(size_t size, const std::type_info &type_info, std::string name);

public:
    FieldBase(Row *row, size_t size, const std::type_info &type_info, std::string name);

    FieldBase(const FieldBase &other);

    FieldBase &operator=(const FieldBase &other);

    FieldBase(FieldBase &&other);

    FieldBase &operator=(FieldBase &&other);

    bool is_comparable(const FieldBase &other) const;

    void add_to_row(Row *row);

    bool belongs_to_row() const;

    bool belongs_to_row(const Row* row) const;

    bool belongs_to_row(const Row& row) const;

    char *begin() const;

    char *end() const;

    Row *row_of_copy() const;

    void zero();

    bool is_zero() const;

    bool operator==(const FieldBase &other) const;
    
    bool operator!=(const FieldBase &other) const;

    bool operator<(const FieldBase &other) const;

    bool operator>(const FieldBase &other) const;

    bool operator<=(const FieldBase &other) const;

    bool operator>=(const FieldBase &other) const;

    defs::hash_t hash() const;

    virtual std::string to_string() const = 0;

    virtual void h5_write_attrs(hid_t parent_handle) {

    }

    virtual void save(hdf5::NdDistListWriter &h5list, const size_t &iitem) const {
        h5list.write_h5item_bytes(iitem, begin());
    }

    virtual void save(hdf5::GroupWriter &gw, size_t irank= 0ul) const {
        gw.save(m_name, begin(), {m_size}, {"raw_data"}, irank);
    }

    virtual void load(hdf5::NdDistListReader &h5list, const size_t &iitem) {
        h5list.read_h5item_bytes(iitem, begin());
    }

    virtual void load(hdf5::GroupReader &gr) {
        gr.load(m_name, begin(), {m_size});
    }

    virtual defs::inds h5_shape() const {
        return {};
    }

    virtual std::vector<std::string> h5_dim_names() const {
        return {};
    }

    virtual hid_t h5_type() const {
        return 0;
    }

private:

    int cmp(const FieldBase &other) const;
};

static std::ostream &operator<<(std::ostream &os, const FieldBase &v) {
    os << v.to_string();
    return os;
}

#endif //M7_FIELDBASE_H
