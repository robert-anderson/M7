//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include <src/core/util/utils.h>
#include <src/core/field/Fields.h>
#include "src/defs.h"
#include "src/core/field/TableField.h"
#include "BufferWindow.h"
#include "src/core/field/Flag.h"

struct TableX {
    std::vector<const TableField *> m_fields;
    BufferWindow m_bw;
    size_t m_row_size;
    size_t m_row_dsize;
    size_t m_tight_row_size = 0ul;
    char *m_data;
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;

    bool is_full() const {
        return m_hwm==m_nrow;
    }

    size_t push_back(size_t nrow=1);

    defs::data_t* ptr();

    char *begin();

    char *begin(const size_t &irow);

    size_t add_field(const TableField *field);

    void move(BufferWindow new_bw);

    void clear();

    void clear_row(const size_t &irow);

    size_t bw_dsize() const;

    std::string field_details(size_t width=30) const;

    void print_field_details(size_t width=30) const;

    void print_contents(const defs::inds* ordering=nullptr) const {
        if (ordering && ordering->size()<m_hwm)
            throw std::runtime_error("ordering array is too short.");
        for (size_t iirow=0ul; iirow<m_hwm; ++iirow){
            auto irow = ordering ? (*ordering)[iirow] : iirow;
            for (auto field: m_fields){
                std::cout << irow << ". " << field->to_string(irow)+" ";
            }
            std::cout << "\n";
        }
    }

};


#endif //M7_TABLE_H
