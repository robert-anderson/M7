//
// Created by rja on 05/03/2020.
//

#include "TableNew.h"

std::vector<std::pair<size_t, size_t>> TableNew::field_offsets() {
    std::vector<std::pair<size_t, size_t>> result{};
    for (auto &field: m_fields) result.push_back(field->m_offset);
    return result;
}

defs::inds TableNew::field_nelements() {
    defs::inds result{};
    for (auto field: m_fields) result.push_back(field->m_nelement);
    return result;
}

void TableNew::print(size_t irow) {
    std::cout << irow << " |";
    for (auto field: m_fields) {
        std::cout << field->to_string(irow) << "  ";
    }
    std::cout << std::endl;
}
