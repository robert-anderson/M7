//
// Created by rja on 07/10/2020.
//

#ifndef M7_PRODUCTENUMERATOR_H
#define M7_PRODUCTENUMERATOR_H

#include <src/core/util/defs.h>
#include <src/core/nd/NdArrayFormat.h>
#include "Enumerator.h"

template<size_t n>
class ProductEnumerator : public Enumerator<defs::inds> {
    const NdArrayFormat<n> m_format;
    size_t iflat = 0ul;
public:
    explicit ProductEnumerator(NdArrayFormat<n>&& format, Enumerator *subsequent = nullptr):
            Enumerator<defs::inds>(subsequent), m_format(std::move(format)) {}

    template<typename ...Args>
    ProductEnumerator(Args... args): ProductEnumerator(NdArrayFormat<n>(args...)){}

    virtual bool next_element(defs::inds &result){
        m_format.decode_flat(iflat, result);
        return iflat++<m_format.nelement();
    }
};


#endif //M7_PRODUCTENUMERATOR_H
