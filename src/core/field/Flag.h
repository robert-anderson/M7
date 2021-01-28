//
// Created by rja on 04/11/2020.
//

#ifndef M7_FLAG_H
#define M7_FLAG_H

#include <cstddef>
#include "src/core/field/BitsetSpecifier.h"
#include "src/core/field/Field.h"

struct Table;

struct FlagBase;

/*
 * FlagField inherits from both a bitset NdField, and a FlagSet.
 * FlagSet is extended with named fields, each of which is an
 * NdFlag, but the FlagSet keeps a pointer to the un-templated base
 * (just as Table keeps pointers to ColumnBase, but in practice the
 * member access is done in terms of the templated NdColumn (actually
 * the NdField, but there are no such complications in the Flag case)).
 */

struct FlagSet {
    typedef NdField<BitsetSpecifier, 0ul> field_t;
    field_t *m_bitset_field = nullptr;
    std::vector<const FlagBase *> m_flags;
    mutable FlagSet *m_last_copied = nullptr;

    FlagSet() {}

    /*
     * The FlagSet may be copied independently of, or as a consequence of the
     * copying of an owning Table. In the former case:
     * 1. using the same idiom as in the Table case, the copy ctor sets the
     *    m_last_copied member of the new object to this ptr, so that copied
     *    Flags (through the FlagBase), can update their pointer to FlagSet.
     * and additionally in the latter case:
     * 2. the m_bitset_field ptr must be set using offsets (see the handling of
     *    MappedTable::m_key_field for analogy)
     * Which of these cases we are in depends on whether m_bitset_field is NULL
     */
    FlagSet(const FlagSet &other);

    size_t add_flag(const FlagBase *flag);

    size_t nbit() const;
};

struct FlagBase {
    FlagSet *m_flagset;
    const size_t m_nelement;
    const size_t m_offset;
    const std::string m_description;

    FlagBase(FlagSet *flagset, size_t nelement, std::string description);

    /*
     * assume that an owning FlagSet was copied
     */
    FlagBase(const FlagBase &other) :
            m_flagset(other.m_flagset->m_last_copied),
            m_nelement(other.m_nelement),
            m_offset(other.m_offset),
            m_description(other.m_description) {}

    BitsetSpecifier::View::BitView operator()(const size_t &irow, const size_t &ielement);

    const BitsetSpecifier::View::BitView operator()(const size_t &irow, const size_t &ielement) const;
};

template<size_t nind>
struct NdFlag : FlagBase {
    NdFormat<nind> m_format;

    template<typename ...Args>
    NdFlag(FlagSet *flagset, std::string description, Args... shape):
            FlagBase(flagset, NdFormat<nind>(shape...).nelement(), description), m_format(shape...) {}

    template<typename ...Args>
    BitsetSpecifier::View::BitView operator()(const size_t &irow, Args... shape) {
        return FlagBase::operator()(irow, m_format.flatten(shape...));
    }

    template<typename ...Args>
    const BitsetSpecifier::View::BitView operator()(const size_t &irow, Args... shape) const {
        return FlagBase::operator()(irow, m_format.flatten(shape...));
    }
};

template<size_t nind>
using Flags = NdFlag<nind>;
using Flag = Flags<0ul>;

#endif //M7_FLAG_H
