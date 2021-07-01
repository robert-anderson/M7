//
// Created by rja on 03/03/2021.
//

#ifndef M7_MEVTABLE_H
#define M7_MEVTABLE_H

#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
/**
 * Multidimensional Expectation Values
 */

template <typename T>
struct MevRow : public Row {
    fields::FermionMevInds m_inds;
    fields::Numbers<T, 1> m_values;

    fields::FermionMevInds &key_field() {
        return m_inds;
    };

    MevRow(size_t nann, size_t ncre, size_t nvalue):
    m_inds(this, nann, ncre, "spin_orbs"), m_values(this, {nvalue}, "values"){}


    MevRow(size_t nop, size_t nvalue): MevRow(nop, nop, nvalue){}
};


#endif //M7_MEVTABLE_H
