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
    const size_t m_element_size;
    const std::type_info &m_type_info;
    const bool m_is_key;
    /*
     * set after instantiation when incorporated into
     * NdField and Row
     */
    RowZ *m_row = nullptr;
    size_t m_nelement = ~0ul;
    size_t m_size = ~0ul;
    size_t m_max_view_offset = ~0ul;
    size_t m_row_offset = ~0ul;

    /*
     * the byte-offset from the beginning of the field corresponding to the current view
     */
    mutable size_t m_view_offset = ~0ul;

    FieldBaseZ(size_t element_size, const std::type_info &type_info, bool is_key);

    FieldBaseZ(const FieldBaseZ &other);

    bool is_added_to_row() const;

    template<size_t nind, typename ...Args> friend struct NdMultiFieldY;
private:

    bool oob() const {
        return m_view_offset>=m_max_view_offset;
    }

    void restart() const;

    void step() const;

    void jump(const size_t& iflat) const;

    bool try_restart() const;

    bool try_step() const;

    bool try_jump(const size_t& iflat) const;

public:
    FieldBaseZ &copy(const FieldBaseZ &other);

    FieldBaseZ &copy_all(const FieldBaseZ &other);

    FieldBaseZ &operator=(const FieldBaseZ &other);

    char *begin() const;

    char *raw_view() const;

    void zero();

    void zero_all();

    bool equals(const FieldBaseZ &other) const;

    bool is_same_type_as(const FieldBaseZ &other) const;

    /**
     * Hashes whole field in row, not the currently viewed element
     * @return
     */
    defs::hash_t hash() const;

    virtual std::string to_string() const = 0;

    std::string to_string_all() const;
};

#endif //M7_FIELDBASEZ_H
