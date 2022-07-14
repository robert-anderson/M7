//
// Created by rja on 14/07/22.
//

#ifndef M7_SPARSE_INVERSE_H
#define M7_SPARSE_INVERSE_H

#include "Fixed.h"
namespace sparse {

    /**
     * more space-efficient representations of the inverse mappings of sparse structures than are available from the
     * dense::Matrix
     */
    namespace inverse {

        template<typename T>
        struct Generic {
            typedef T elem_t;
            static_assert(std::is_base_of<Element, T>::value, "template arg must be derived from Element");
            typedef std::pair<uint_t, uint_t> key_t;
            typedef std::pair<key_t, T> pair_t;
        private:
            const T m_not_found_entry;
            std::map<key_t, T> m_conns;
        public:
            Generic(const dynamic::Generic<T>& src) {
                for (uint_t irow = 0ul; irow < src.nrow(); ++irow) {
                    for (auto& it: src[irow]) {
                        const auto elem = static_cast<const Element&>(it);
                        m_conns.insert({{irow, elem.m_i}, it});
                    }
                }
            }

            Generic(const fixed::Generic<T>& src) {
                for (uint_t irow = 0ul; irow < src.m_nrow; ++irow) {
                    for (auto it = src.cbegin(irow); it != src.cend(irow); ++it) {
                        const auto elem = static_cast<const Element&>(*it);
                        m_conns.insert({{irow, elem.m_i}, *it});
                    }
                }
            }

            const T& entry(uint_t i, uint_t j) const {
                auto it = m_conns.find({i, j});
                return it == m_conns.cend() ? m_not_found_entry : it->second;
            }

            bool exists(uint_t i, uint_t j) const {
                return m_conns.find({i, j}) != m_conns.cend();
            }
        };

        typedef Generic<Element> Network;

        template<typename T>
        struct Matrix : Generic<MatrixElement<T>> {
            typedef MatrixElement<T> elem_t;

            Matrix(const dynamic::Matrix<T>& src) : Generic<elem_t>(src) {}

            Matrix(const fixed::Matrix<T>& src) : Generic<elem_t>(src) {}

            T get(uint_t i, uint_t j) const {
                return Generic<elem_t>::entry(i, j).m_v;
            }
        };
    } // inverse
} // sparse

#endif //M7_SPARSE_INVERSE_H
