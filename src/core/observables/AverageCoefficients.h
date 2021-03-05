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
    buffered::Vector<defs::mev_ind_t> m_working_inds;
    defs::wf_t m_ref_coeff = 0.0;

    AverageCoefficients(std::string name, defs::inds ninds, size_t nvalue, size_t nbucket = 100) :
            BufferedTable<MevRow<defs::wf_t>, true>(name, {{ninds, nvalue}, nbucket}),
            m_working_inds(m_row.nelement()) {}


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
        size_t i = 0ul;
        for (size_t j = 0ul; j < key.ncre(); ++j) m_working_inds(i++) = key.cre(j);
        for (size_t j = 0ul; j < key.nann(); ++j) m_working_inds(i++) = key.ann(j);
    }
};


#endif //M7_AVERAGECOEFFICIENTS_H
