//
// Created by RJA on 21/09/2020.
//

#ifndef M7_NDARRAY_H
#define M7_NDARRAY_H

#include <string>
#include <vector>
#include "NdSpecifier.h"
#include "src/core/util/defs.h"

/*
 * efficient header-only multidimensional array implementation for array-like objects
 * with a compile time-determined number of dimensions, but a constant, run
 * time-determined shape.
 */

namespace nd_array {

    template<typename T>
    struct Accessor {
        T &m_data;

        Accessor(T &data) : m_data(data) {}

        std::string to_string() const {
            return std::to_string(m_data);
        }

        operator T &() { return m_data; }

        operator const T &() const { return m_data; }

        T &operator=(const T &v) {
            m_data = v;
            return *this;
        }
    };

    template<typename T, size_t nind>
    struct Selector {
        typedef Accessor<T> accessor_t;
        typedef Accessor<const T> const_accessor_t;
        const NdSpecifier<Selector, nind> &m_format;
        std::vector <T> &m_data;

        Selector(
                const NdSpecifier<Selector, nind> &format,
                std::vector <T> &data) : m_format(format), m_data(data) {}

        accessor_t operator()(const size_t &flat) {
            ASSERT(flat < m_format.nelement());
            return accessor_t(m_data[flat]);
        }

        const_accessor_t operator()(const size_t &flat) const {
            ASSERT(flat < m_format.nelement());
            return const_accessor_t(m_data[flat]);
        }
    };


    template<typename T, size_t nind>
    class Specifier : public NdSpecifier<Selector<T, nind>, nind> {
        std::vector <T> m_data;
        typedef NdSpecifier<Selector<T, nind>, nind> base_t;
    public:

        using base_t::nelement;

        template<typename ...Args>
        Specifier(const size_t &first, Args ...shape):
                base_t(Selector<T, nind>(*this, m_data), first, shape...),
                m_data(nelement(), 0) {}
    };

}



#endif //M7_NDARRAY_H
