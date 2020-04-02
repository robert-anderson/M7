//
// Created by Robert John Anderson on 2020-03-26.
//

#include "Flag.h"
#include "FlagField.h"


Flag::Flag(FlagField* field, size_t nelement):
m_field(field), m_nelement(nelement), m_offset(field->add_flag(this)) {}

FlagElement Flag::element(const size_t &irow, const size_t &isegment, const size_t &ielement) {
    assert(ielement<m_nelement);
    return FlagElement(m_field->element(irow, isegment), ielement);
}
