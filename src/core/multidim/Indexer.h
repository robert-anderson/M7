//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_INDEXER_H
#define M7_INDEXER_H

#include <cstddef>
#include <array>
#include <iostream>
#include "src/utils.h"

template<size_t nind>
class Indexer {
protected:
    std::array<size_t, nind> m_shape{};
    std::array<size_t, nind> m_strides{};
public:

    template<typename ...Args>
    Indexer(const size_t &first, Args ...shape) {
        static_assert(sizeof...(shape) + 1 == nind, "Invalid number of arguments.");
        set_shape(first, shape...);
        set_strides();
    }

    Indexer(const std::array<size_t, nind> &shape) {
        m_shape = shape;
        set_strides();
    }

    const std::array<size_t, nind> &shape() const {
        return m_shape;
    }

    const std::array<size_t, nind> &strides() const {
        return m_strides;
    }

    size_t nelement() const {
        return m_strides.front() * m_shape.front();
    }

private:
    template<typename T>
    void set_shape(const T &extent) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape.fill(extent);
    }

    template<typename T, typename ...Args>
    void set_shape(const T &first, Args ...args) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind - sizeof...(args) - 1] = first;
        set_shape(args...);
    }

    template<typename T>
    void set_shape(const T &first, T &second) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind - 2] = first;
        m_shape[nind - 1] = second;
    }

    void set_strides() {
        m_strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            m_strides[nind - i] = m_strides[nind - i + 1] * m_shape[nind - i + 1];
        }
    }

    template<size_t nind_get>
    size_t partial_get() const {
        return 0;
    }

    template<size_t nind_get>
    size_t partial_get(const size_t &first) const {
        return m_strides[nind_get - 1] * first;
    }

    template<size_t nind_get, typename ...Args>
    size_t partial_get(const size_t &first, Args ...trailing) const {
        ASSERT(first < m_shape[nind - sizeof...(trailing) - 1]);
        constexpr size_t iind = nind_get - sizeof...(trailing) - 1;
        if (iind < nind_get) {
            return m_strides[iind] * first + partial_get<nind_get>(trailing...);
        } else {
            return m_strides[iind] * first;
        }
    }


public:
    template<typename ...Args>
    size_t get_sub(Args...inds) const {
        static_assert(sizeof...(inds) <= nind, "Invalid number of arguments.");
        return partial_get<sizeof...(inds)>(inds...);
    }

    template<typename ...Args>
    size_t get(Args ...inds) const {
        static_assert(sizeof...(inds) == nind, "Invalid number of arguments.");
        return get_sub(inds...);
    }
};

#endif //M7_INDEXER_H
