//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDBASEZ_H
#define M7_FIELDBASEZ_H

#include <src/defs.h>
#include <cstring>
#include <src/core/hash/Hashing.h>
#include "RowZ.h"

struct FieldBaseZ {
    RowZ *m_row;
    const size_t m_item_size;
    const std::type_info &m_type_info;
    const size_t m_nitem;
    const size_t m_size;
    size_t m_row_offset = ~0ul;

    std::vector<char> m_null_string;

    FieldBaseZ(RowZ* row, size_t item_size, size_t nitem, const std::type_info &type_info);

    FieldBaseZ(const FieldBaseZ &other);

    FieldBaseZ& operator=(const FieldBaseZ& other);

    bool is_comparable(const FieldBaseZ& other) const;

    void add_to_row(RowZ* row);

    bool is_added_to_row() const;

    char *begin() const;

    char *begin(const size_t& ielement) const;

    char *end() const;

    char *end(const size_t& ielement) const;

    RowZ* row_of_copy() const;

    void zero();

    bool is_zero() const;

    bool operator==(const FieldBaseZ &other) const;

    /**
     * Hashes whole field in row, not the currently viewed element
     * @return
     */
    defs::hash_t hash() const;

    virtual std::string to_string_element(const size_t& iitem) const = 0;

    std::string to_string() const;
};


#endif //M7_FIELDBASEZ_H
