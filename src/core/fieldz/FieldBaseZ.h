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

    std::vector<char> m_null_item_string;

    FieldBaseZ(RowZ* row, size_t item_size, size_t nitem, const std::type_info &type_info);

    FieldBaseZ(const FieldBaseZ &other);

    FieldBaseZ& operator=(const FieldBaseZ& other){
        if (&other==this) return *this;
        MPI_ASSERT(is_comparable(other),
                   "can't copy from a field which is either incompatible or has a different selection length")
        std::memcpy(begin(), other.begin(), m_size);
        return *this;
    }

    bool is_comparable(const FieldBaseZ& other) const {
        return m_item_size==other.m_nitem &&
               m_type_info==other.m_type_info && m_size==other.m_size;
    }

    bool is_added_to_row() const;

    //template<size_t nind, typename ...Args> friend struct NdMultiFieldZ;

    char *begin() const;

    char *begin(const size_t& ielement) const;

    char *end() const;

    char *end(const size_t& ielement) const;

    void zero();

    bool is_zero() const;

    bool operator==(const FieldBaseZ &other) const {
        MPI_ASSERT(is_comparable(other),
                   "can't copy from a field which is either incompatible or has a different selection length")
        return std::memcmp(begin(), other.begin(), m_size) == 0;
    }

    /**
     * Hashes whole field in row, not the currently viewed element
     * @return
     */
    defs::hash_t hash() const;

    virtual std::string to_string_element(const size_t& iitem) const = 0;

    std::string to_string() const;
};


#endif //M7_FIELDBASEZ_H
