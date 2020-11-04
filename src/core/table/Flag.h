//
// Created by rja on 04/11/2020.
//

#ifndef M7_FLAG_H
#define M7_FLAG_H

#include <cstddef>
#include "src/core/field/BitsetSpecifier.h"
#include "src/core/nd/NdFormat.h"

struct TableX;

struct TableFlag {
    TableX* m_table;
    const size_t m_nelement;
    const size_t m_offset;
    TableFlag(TableX* table, size_t nelement);

    BitsetSpecifier::View::BitView operator()(const size_t& irow, const size_t& ielement);
};

template <size_t nind>
struct NdFlag : TableFlag {
    NdFormat<nind> m_format;

    template<typename ...Args>
    NdFlag(TableX* table, Args... shape):
    TableFlag(table, NdFormat<nind>(shape...).nelement()), m_format(shape...){}

    template<typename ...Args>
    BitsetSpecifier::View::BitView operator()(const size_t& irow, Args... shape){
        return TableFlag::operator()(irow, m_format.flatten(shape...));
    }
};


#endif //M7_FLAG_H
