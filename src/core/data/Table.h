//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <src/core/util/utils.h>
#include "src/core/util/defs.h"
#include "NdFieldBase.h"
#include "BufferWindow.h"

struct TableX {
    BufferWindow m_bw;
    size_t m_row_size;
    size_t m_row_dsize;
    size_t m_tight_row_size = 0ul;
    std::vector<const NdFieldBaseX *> m_fields;
    char *m_data;
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;


    //TableX(BufferWindow bw): m_bw(bw){}

    char *begin() {
        return (char *) m_bw.m_ptr;
    }

    char *begin(const size_t &irow) {
        return begin() + irow * m_row_size;
    }

    size_t add_field(const NdFieldBaseX *field) {
        // returns the offset in bytes for the field being added
        auto offset = 0ul;
        if(!m_fields.empty()){
            offset = m_fields.back()->m_offset+m_fields.back()->m_size;
            if (!m_fields.back()->is_same_type_as(*field)){
                // go to next whole dataword
                offset = integer_utils::divceil(offset, defs::nbyte_data)*defs::nbyte_data;
            }
        }

        m_tight_row_size = offset+field->m_size;
        m_row_dsize = integer_utils::divceil(m_tight_row_size, defs::nbyte_data);
        m_row_size = m_row_dsize*defs::nbyte_data;

        m_fields.push_back(field);
        return offset;
    }

    void move(BufferWindow new_bw) {
        if (m_bw) std::memmove(m_bw.m_ptr, new_bw.m_ptr, sizeof(defs::data_t) * std::min(m_bw.m_dsize, new_bw.m_dsize));
        m_bw = new_bw;
        if (!m_row_size) return;
        m_nrow = (sizeof(defs::data_t) * m_bw.m_dsize) / m_row_size;
    }

    void clear() {
        std::memset((char *) (m_bw.m_ptr), 0, m_row_size * m_hwm);
        m_hwm = 0ul;
    }

    void clear_row(const size_t &irow) {
        std::memset(begin(irow), 0, m_tight_row_size);
    }

    std::string field_details(size_t width=30) const {
        std::string res;
        for (size_t i = 0ul; i < m_fields.size(); ++i) {
            std::string desc;
            if (!m_fields[i]->m_description.empty()) desc = " (\""+m_fields[i]->m_description+"\")";
            res += "\nField " + std::to_string(i) + desc + ":\n";
            for (auto pair:m_fields[i]->m_details)
                res += "\t" + utils::padded_string(pair.first, width) + ": " + utils::padded_string(pair.second, width)+"\n";
        }
        return res;
    }

    void print_field_details(size_t width=30) const {
        std::cout << field_details(width) << std::endl;
    }

};


#endif //M7_TABLE_H
