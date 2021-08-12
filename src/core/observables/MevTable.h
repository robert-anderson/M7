//
// Created by rja on 03/03/2021.
//

#ifndef M7_MEVTABLE_H
#define M7_MEVTABLE_H

#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"

/**
 * Multidimensional Averaged Estimators
 */

struct MaeInds : fields::Numbers<T, 1>


template <typename T>
struct MaeRow : public Row {
    fields::FermionMevInds m_inds;
    fields::Numbers<T, 1> m_values;

    fields::FermionMevInds &key_field() {
        return m_inds;
    };

    MaeRow(size_t nann, size_t ncre, size_t nvalue):
    m_inds(this, nann, ncre, "spin_orbs"), m_values(this, {nvalue}, "values"){}


    MaeRow(size_t nop, size_t nvalue): MaeRow(nop, nop, nvalue){}
};


#endif //M7_MEVTABLE_H
