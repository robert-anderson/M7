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
class MevRow : public Row {
    defs::inds get_offsets(const defs::inds& ninds){
        defs::inds tmp;
        tmp.reserve(ninds.size());
        tmp.emplace_back(0);
        for (size_t i=1ul; i < ninds.size(); ++i) tmp.push_back(tmp[i - 1] + ninds[i]);
        return tmp;
    }

public:

    defs::inds m_ninds;
    defs::inds m_offsets;
    fields::Vector<defs::mev_ind_t> m_inds;
    fields::Vector<T> m_values;

    fields::Vector<defs::mev_ind_t> &key_field() {
        return m_inds;
    };

    MevRow(defs::inds ninds, size_t nvalue):
    m_ninds(ninds), m_offsets(get_offsets(ninds)),
    m_inds(this, nelement()), m_values(this, nvalue){}

    size_t nelement() const {
        return m_offsets.back()+m_ninds.back();
    }

    void set(const size_t& iind, const defs::inds& inds){
        ASSERT(iind<m_ninds.size());
        ASSERT(inds.size()==m_ninds[iind]);
        const auto offset = m_offsets[iind];
        // implicitly convert to storage type defs::mev_ind_t
        for (size_t i = 0ul; i<inds.size(); ++i) {
            m_values(offset+i) = inds[i];
        }
    }

    void get(const size_t& iind, defs::inds& inds) const{
        ASSERT(iind<m_ninds.size());
        ASSERT(inds.size()==m_ninds[iind]);
        const auto offset = m_offsets[iind];
        // implicitly convert to storage type defs::mev_ind_t
        for (size_t i = 0ul; i<inds.size(); ++i) inds[i] = m_values(offset+i);
    }
};


#endif //M7_MEVTABLE_H
