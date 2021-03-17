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

struct FieldBase {
    Row *m_row;
    const size_t m_item_size;
    const std::type_info &m_type_info;
    const size_t m_nitem;
    const size_t m_size;
    size_t m_row_offset = ~0ul;

    std::vector<char> m_null_string;
    hid_t m_h5type;

    FieldBase(Row* row, size_t item_size, size_t nitem, const std::type_info &type_info, hid_t h5type);

    FieldBase(const FieldBase &other);

    FieldBase& operator=(const FieldBase& other);

    bool is_comparable(const FieldBase& other) const;

    void add_to_row(Row* row);

    bool is_added_to_row() const;

    char *begin() const;

    char *begin(const size_t& ielement) const;

    char *end() const;

    char *end(const size_t& ielement) const;

    Row* row_of_copy() const;

    void zero();

    bool is_zero() const;

    bool operator==(const FieldBase &other) const;

    bool operator!=(const FieldBase &other) const;

    bool operator<(const FieldBase &other) const;

    bool operator>(const FieldBase &other) const;

    bool operator>=(const FieldBase &other) const;

    bool operator<=(const FieldBase &other) const;

    /**
     * Hashes whole field in row, not the currently viewed element
     * @return
     */
    defs::hash_t hash() const;

    virtual std::string to_string_element(const size_t& iitem) const = 0;

    std::string to_string() const;

    virtual void h5_write(hdf5::NdListWriterBase& h5list, const size_t& iitem){
        h5list.write_h5item_bytes(iitem, begin());
    }

    virtual void h5_read(hdf5::NdListReaderBase& h5list, const size_t& iitem){
        h5list.read_h5item_bytes(iitem, begin());
    }

private:
    int cmp(const FieldBase &other) const;
};


#endif //M7_FIELDBASE_H
