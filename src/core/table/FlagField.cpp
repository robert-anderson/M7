//
// Created by Robert John Anderson on 2020-03-26.
//

#include "FlagField.h"
#include "Table.h"


defs::pair FlagField::add_flag(Flag *flag) {
    defs::pair offset;
    if (m_flags.empty()) offset={0,0};
    else {
        offset.first = m_flags.back()->m_offset.first;
        offset.second = m_flags.back()->m_offset.second+m_flags.back()->m_nelement;
    }
    m_flags.push_back(flag);
    increment_nbit(flag->m_nelement);
    m_table->update_last_field();
    return rectify_offset(offset);
}