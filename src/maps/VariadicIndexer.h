//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_VARIADICINDEXER_H
#define M7_VARIADICINDEXER_H

#include <cstddef>
#include <array>
#include <iostream>
#include "../utils.h"

template <size_t nind>
class VariadicIndexer {
    std::array<size_t, nind> m_shape{};
    std::array<size_t, nind> m_strides{};
public:
    template<typename ...Args>
    VariadicIndexer(Args ...args) {
        static_assert(sizeof...(args)==1 || sizeof...(args)==nind, 
				"Invalid number of arguments.");
        set_shape(args...);
        m_strides[nind-1] = 1ul;
        for (auto i=2ul; i<=nind; i++){
            m_strides[nind-i] = m_strides[nind-i+1]*m_shape[nind-i];
        }
        utils::print(m_strides);
    }

private:
    template<typename T>
    std::array<size_t, nind> set_shape(const T &extent){
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape.fill(extent);
    }
    template<typename T, typename ...Args>
    std::array<size_t, nind> set_shape(const T &first, Args ...args){
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind-sizeof...(args)-1] = first;
        set_shape(args...);
    }
    template<typename T>
    std::array<size_t, nind> set_shape(const T &first, T &second){
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind-2] = first;
        m_shape[nind-1] = second;
    }

    template<typename T>
    size_t partial_get(const T &first){
        static_assert(std::is_integral<T>::value, "Indexing requires an integral value.");
        assert(first<m_shape.back());
        return first;
    }

    template<typename T, typename ...Args>
    size_t partial_get(const T &first, Args ...args){
        static_assert(std::is_integral<T>::value, "Indexing requires an integral value.");
        assert(first<m_shape[nind-sizeof...(args)-1]);
        return m_strides[nind-sizeof...(args)-1]*first+partial_get(args...);
    }

public:
    template<typename ...Args>
    size_t get(Args ...args){
        static_assert(sizeof...(args)==nind, "Invalid number of arguments.");
        return partial_get(args...);
    }
};


#endif //M7_VARIADICINDEXER_H
