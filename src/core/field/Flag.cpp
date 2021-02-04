//
// Created by rja on 04/11/2020.
//

#include "Flag.h"

#include <utility>
#include "src/core/table/Table.h"


size_t FlagSet::add_flag(const FlagBase *flag) {
    auto offset = m_flags.empty() ? 0ul : m_flags.back()->m_offset+m_flags.back()->m_nelement;
    m_flags.push_back(flag);
    return offset;
}

size_t FlagSet::nbit() const {
    return m_flags.empty() ? 0ul : m_flags.back()->m_offset+m_flags.back()->m_nelement;
}

FlagSet::FlagSet(const FlagSet &other) {
    if (other.m_bitset_field){
        // Owning Table was copied
        auto src_table = other.m_bitset_field->m_column.m_table;
        auto cpy_table = src_table->m_last_copied;
        m_bitset_field = (field_t*)((const char*)cpy_table + other.m_bitset_field->symbol_offset(src_table));
    }
    other.m_last_copied = this;
}

FlagBase::FlagBase(FlagSet *flagset, size_t nelement, std::string description) :
        m_flagset(flagset),
        m_nelement(nelement),
        m_offset(flagset->add_flag(this)),
        m_description(std::move(description)){}

BitsetSpecifier::View::BitView FlagBase::operator()(const size_t &irow, const size_t &ielement) {
    ASSERT(ielement<m_nelement);
    return (*m_flagset->m_bitset_field)(irow) [m_offset+ielement];
}

const BitsetSpecifier::View::BitView FlagBase::operator()(const size_t &irow, const size_t &ielement) const {
    ASSERT(ielement<m_nelement);
    return (*m_flagset->m_bitset_field)(irow) [m_offset+ielement];
}
