//
// Created by rja on 05/03/2020.
//

#include "FieldSet.h"

std::vector<std::pair<size_t, size_t>> FieldSet::field_offsets() {
    std::vector<std::pair<size_t, size_t>> result{};
    for (auto &field: m_fields) result.push_back(field->m_offset);
    return result;
}

defs::inds FieldSet::field_nelements() {
    defs::inds result{};
    for (auto field: m_fields) result.push_back(field->m_nelement);
    return result;
}

void FieldSet::print(size_t irow) {
    std::cout << irow << " |";
    for (auto field: m_fields) {
        for (size_t ientry = 0ul; ientry < field->m_nelement; ++ientry) {
            std::cout << field->to_string(irow, ientry);
            std::cout << "   ";
        }
    }
    std::cout << std::endl;
}
