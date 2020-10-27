//
// Created by RJA on 27/10/2020.
//

#ifndef M7_NDFIELDBASE_H
#define M7_NDFIELDBASE_H

#include "Field.h"

struct TableX;

struct NdFieldBaseX : FieldBaseX {
    TableX* m_table;
    const std::string m_description;
    const size_t m_nelement;
    const size_t m_size;
    const size_t m_offset;

    NdFieldBaseX(TableX* table, size_t nelement, size_t element_size,
        const std::type_info &type_info, std::string description);

    char* begin(const size_t &irow) const;

    char* raw_ptr(const size_t &irow, const size_t &ielement) const;

    std::string to_string(size_t irow) const {
        std::string res;
        for (size_t ielement=0ul; ielement<m_nelement; ++ielement)
            res+=element_string(raw_ptr(irow, ielement))+" ";
        return res;
    }

    std::map<std::string, std::string> details() const override;

};


#endif //M7_NDFIELDBASE_H
