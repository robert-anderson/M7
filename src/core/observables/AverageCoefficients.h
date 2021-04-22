//
// Created by rja on 03/03/2021.
//

#ifndef M7_AVERAGECOEFFICIENTS_H
#define M7_AVERAGECOEFFICIENTS_H

#include <src/core/basis/Connections.h>
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/BufferedFields.h"
#include "MevTable.h"

struct AverageCoefficients : BufferedTable<MevRow<defs::wf_t>, true> {
    buffered::FermionMevInds m_working_inds;
    defs::wf_t m_ref_coeff = 0.0;

    AverageCoefficients(std::string name, size_t nann, size_t ncre, size_t nvalue, size_t nbucket = 100) :
            BufferedTable<MevRow<defs::wf_t>, true>(name, {{nann, ncre, nvalue}, nbucket}),
            m_working_inds({nann, ncre}) {}


    LookupResult operator[](const conn::Basic<0> &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::operator[](m_working_inds);
    }

    size_t insert(const conn::Basic<0> &key) {
        set_working_inds(key);
        return MappedTable<MevRow<defs::wf_t>>::insert(m_working_inds);
    }

private:
    void set_working_inds(const conn::Basic<0> &key) {
        m_working_inds = {key.ann(), key.cre()};
    }
};


#endif //M7_AVERAGECOEFFICIENTS_H
