//
// Created by rja on 07/10/2020.
//

#ifndef M7_PRODUCTENUMERATOR_H
#define M7_PRODUCTENUMERATOR_H

#include <M7_lib/defs.h>
#include <M7_lib/nd/NdFormat.h>

#include "Enumerator.h"

class ProductEnumerator : public Enumerator<defs::inds> {
    const defs::inds m_shape;
    const size_t m_nind;
    const defs::inds m_strides;
    const size_t m_nelement;
    size_t m_iflat = 0ul;

    defs::inds strides(){
        defs::inds res;
        if (!m_nind) return {};
        res.resize(m_nind);
        res.back() = 1ul;
        for (auto i = 2ul; i <= m_nind; i++) {
            res[m_nind - i] = res[m_nind - i + 1] * m_shape[m_nind - i + 1];
        }
        return res;
    }

    void decode_flat(const size_t& iflat, defs::inds& inds) const {
        ASSERT(inds.size()>=m_nind);
        size_t remainder = iflat;
        for (size_t i=0ul; i<m_nind; ++i){
            auto& ind = inds[i];
            ind = remainder/m_strides[i];
            remainder-=ind*m_strides[i];
        }
    }

public:

    explicit ProductEnumerator(defs::inds&& shape, Enumerator *subsequent = nullptr):
            Enumerator<defs::inds>(subsequent), m_shape(std::move(shape)),
            m_nind(m_shape.size()), m_strides(strides()),
            m_nelement(m_shape.size()?m_shape.front()*m_strides.front():1){}

    ProductEnumerator(size_t nind, size_t extent, Enumerator *subsequent = nullptr):
    ProductEnumerator(defs::inds(nind, extent), subsequent){}

    virtual bool next_element(defs::inds &result){
        decode_flat(m_iflat, result);
        auto allfound = m_iflat++==m_nelement;
        if (allfound) m_iflat = 0ul;
        return !allfound;
    }
};


#endif //M7_PRODUCTENUMERATOR_H
