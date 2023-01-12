//
// Created by Robert J. Anderson on 14/06/2020.
//

#ifndef M7_SPARSE_H
#define M7_SPARSE_H

#include <vector>
#include <forward_list>

#include "M7_lib/parallel/MPIAssert.h"
#include "M7_lib/io/Logging.h"
#include <set>
#include "M7_lib/util/Convert.h"

namespace sparse {

    /**
     * In a sparse structure, rows are made up of many elements.
     * each row corresponds to a single basis vector of a matrix, or node of a network
     * A network element carries no more than an integer identifying the connected node index.
     */
    struct Element {
        /**
         * sparsely-connected index
         */
        uint_t m_i;
        Element(uint_t i=~0ul): m_i(i){}

        virtual ~Element() = default;

        virtual str_t to_string() const {
            return convert::to_string(m_i);
        }

        operator uint_t () const {
            return m_i;
        }

        template<typename T>
        static const Element& cast(const T& elem) {
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            return elem;
        }

        template<typename T>
        static Element& cast(T& elem) {
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            return elem;
        }
    };

    template<typename T>
    struct MatrixElement : Element {
        /**
         * value associated with the sparsely-connected element
         */
        T m_v;

        MatrixElement(uint_t i=~0ul, T v={}): Element(i), m_v(v){}

        str_t to_string() const override {
            return logging::format("({} -> {})", m_i, convert::to_string(m_v));
        }

        operator str_t () const {
            return to_string();
        }

        bool operator==(const MatrixElement<T>& other) const {
            return m_i==other.m_i && m_v==other.m_v;
        }

        bool operator!=(const MatrixElement<T>& other) const {
            return !(*this==other);
        }
    };
}

#endif //M7_SPARSE_H