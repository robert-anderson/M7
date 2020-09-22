//
// Created by RJA on 21/09/2020.
//

#ifndef M7_NDINDICES_H
#define M7_NDINDICES_H

#include "NdSpecifier.h"
#include <cstddef>

/*
 * trivial case, does not map to data buffer, just converts multidimensional indices
 * to flat index.
 */

namespace nd_indices {

    template<size_t nind>
    struct Selector {
        typedef size_t accessor_t;
        typedef size_t const_accessor_t;

        accessor_t operator()(const size_t &flat) {
            return flat;
        }

        const_accessor_t operator()(const size_t &flat) const {
            return flat;
        }
    };

    template<size_t nind>
    class Specifier : public NdSpecifier<Selector<nind>, nind> {
        typedef NdSpecifier<Selector<nind>, nind> base_t;
    public:

        using base_t::nelement;

        template<typename ...Args>
        Specifier(const size_t &first, Args ...shape):
                base_t(Selector<nind>(), first, shape...) {}
    };


}

#endif //M7_NDINDICES_H
