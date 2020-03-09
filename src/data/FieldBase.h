//
// Created by Robert John Anderson on 2020-03-09.
//

#ifndef M7_FIELDBASE_H
#define M7_FIELDBASE_H

#include <string>

struct FieldSet;

struct FieldBase {
    FieldSet *m_field_set;
    const size_t m_length;
    /*
     * offset in number of field entries stored from the beginning of the row
     */
    const size_t m_offset = 0ul;

    template<typename ...Args>
    std::string to_string(size_t irow, Args... inds){return "";};

    FieldBase(FieldSet *field_set, size_t length);
};


#endif //M7_FIELDBASE_H
