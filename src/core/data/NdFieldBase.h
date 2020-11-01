//
// Created by RJA on 27/10/2020.
//

#ifndef M7_NDFIELDBASE_H
#define M7_NDFIELDBASE_H

#include "Field.h"
#include "src/core/hash/Hashing.h"

struct TableX;

struct NdFieldBaseX : FieldBaseX {
    TableX *m_table;
    const std::string m_description;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_offset;

//    NdFieldBaseX(TableX *table, size_t nelement, size_t element_size,
//                 const std::type_info &type_info, std::string description);

    NdFieldBaseX(TableX *table, FieldBaseX &&field,
                 size_t nelement, std::string description);

    char *begin(const size_t &irow) const;

    char *raw_ptr(const size_t &irow, const size_t &ielement) const {
        return begin(irow) + ielement * m_element_size;
    }

//    typedef std::pair<const char *, size_t> raw_view_t;


//    struct equals_fn {
//        bool operator()(const raw_view_t& v1, const raw_view_t& v2){
//            ASSERT(v1.second==v2.second);
//            return !std::memcmp(v1.first, v2.first, v1.second);
//        }
//    };

//    struct default_hash_fn {
//        defs::hash_t operator()(const raw_view_t& v){
//            return hashing::fnv_hash(v.first, v.second);
//        }
//    };

    std::string to_string(size_t irow) const {
        std::string res;
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement)
            res += element_string(raw_ptr(irow, ielement)) + " ";
        return res;
    }
};


#endif //M7_NDFIELDBASE_H
