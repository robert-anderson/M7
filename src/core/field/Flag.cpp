//
// Created by rja on 04/11/2020.
//

#include "Flag.h"
#include "src/core/table/Table.h"

TableFlag::TableFlag(TableX *table, size_t nelement) :
        m_table(table),
        m_nelement(nelement),
        m_offset(table->add_flag(this)){}

BitsetSpecifier::View::BitView TableFlag::operator()(const size_t &irow, const size_t &ielement) {
    ASSERT(ielement<m_nelement);
    return (*m_table->m_flag_field)(irow)[m_offset+ielement];
}
