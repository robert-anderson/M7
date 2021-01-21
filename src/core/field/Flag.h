//
// Created by rja on 04/11/2020.
//

#ifndef M7_FLAG_H
#define M7_FLAG_H

#include <cstddef>
#include "src/core/field/BitsetSpecifier.h"
#include "src/core/field/TableField.h"

struct Table;

struct FlagBase;

struct FlagSet {
    NdField<BitsetSpecifier, 0ul> * m_bitset_field;
    FlagSet(NdField<BitsetSpecifier, 0ul>* bitset_field):
    m_bitset_field(bitset_field){}
    std::vector<const FlagBase*> m_flags;
    size_t add_flag(const FlagBase* flag);
    size_t nbit() const;
};

struct FlagBase {
    FlagSet *m_flagset;
    const size_t m_nelement;
    const size_t m_offset;
    const std::string m_description;
    FlagBase(FlagSet *flagset, size_t nelement, std::string description);

    BitsetSpecifier::View::BitView operator()(const size_t& irow, const size_t& ielement);

    const BitsetSpecifier::View::BitView operator() (const size_t& irow, const size_t& ielement) const;
};

template <size_t nind>
struct NdFlag : FlagBase {
    NdFormat<nind> m_format;

    template<typename ...Args>
    NdFlag(FlagSet *flagset, std::string description, Args... shape):
            FlagBase(flagset, NdFormat<nind>(shape...).nelement(), description), m_format(shape...){}

    template<typename ...Args>
    BitsetSpecifier::View::BitView operator()(const size_t& irow, Args... shape){
        return FlagBase::operator()(irow, m_format.flatten(shape...));
    }

    template<typename ...Args>
    const BitsetSpecifier::View::BitView operator()(const size_t& irow, Args... shape) const {
        return FlagBase::operator()(irow, m_format.flatten(shape...));
    }
};

template <size_t nind>
using Flags = NdFlag<nind>;
using Flag = Flags<0ul>;

#endif //M7_FLAG_H
