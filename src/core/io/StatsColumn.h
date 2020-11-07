//
// Created by rja on 07/11/2020.
//

#ifndef M7_STATSCOLUMN_H
#define M7_STATSCOLUMN_H

#include <src/core/enumerator/ProductEnumerator.h>
#include "src/core/nd/NdArray.h"

struct StatsSpecifier;

struct StatsColumnBase {
    const size_t m_nelement;
    const defs::inds m_shape;
    const size_t m_nsubcolumn;
    const std::string m_description;

    StatsColumnBase(StatsSpecifier *spec, size_t nelement, defs::inds shape, size_t nsubcolumn,
                    std::string description);

    ProductEnumerator format_enum() const {
        return ProductEnumerator(defs::inds(m_shape));
    }

    virtual std::string to_string() const = 0;
    virtual void zero() = 0;
};

template<typename T, size_t nind = 0>
struct StatsColumn : NdArray<T, nind>, StatsColumnBase {
    using NdArray<T, nind>::nelement;
    using NdArray<T, nind>::m_data;

    template<typename ...Args>
    StatsColumn(StatsSpecifier *spec, std::string description, Args... shape):
            NdArray<T, nind>(shape...),
            StatsColumnBase(spec, nelement(),
                            defs::inds(
                                    NdArray<T, nind>::m_format.shape().begin(),
                                    NdArray<T, nind>::m_format.shape().end()
                            ),
                            consts::is_complex<T>() ? 2ul : 1ul, description) {}

    std::string to_string() const override {
        std::string res;
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement) {
            res += " " + std::to_string(consts::real(m_data[ielement]));
            if (consts::is_complex<T>())
                res += " " + std::to_string(consts::imag(m_data[ielement]));
        }
        return res;
    }

    void zero() override {
        NdArray<T, nind>::zero();
    }
};


#endif //M7_STATSCOLUMN_H