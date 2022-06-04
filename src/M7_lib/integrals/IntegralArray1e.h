//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY1E_H
#define M7_INTEGRALARRAY1E_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/defs.h>
#include "IntegralStorage.h"

namespace integrals_1e {
    struct IndexerSymNone : IntegralIndexer {
        IndexerSymNone(size_t norb);

        size_t index_only(size_t a, size_t i) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t i) const;
    };

    using namespace integer_utils;

    struct IndexerSymH : IntegralIndexer {
        IndexerSymH(size_t norb);

        size_t index_only(size_t a, size_t i) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t i) const;
    };

    template<typename T>
    struct Array {
        virtual bool set(size_t a, size_t i, T elem) = 0;

        virtual T get(size_t a, size_t i) const = 0;
    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<IntegralIndexer, indexer_t>::value, "invalid template class");
        using Array<T>::set;
        using Array<T>::get;
        indexer_t m_indexer;
        PrivateIntegralStorage<T> m_data;

        IndexedArray(size_t norb) : m_indexer(norb), m_data(static_cast<const IntegralIndexer &>(m_indexer).m_size) {}

        bool set(size_t a, size_t i, T elem) override {
            // any compiler should statically execute this conditional
            if (!consts::is_complex<T>())
                return m_data.set_data(m_indexer.index_only(a, i), elem);
            else {
                const auto pair = m_indexer.index_and_conj(a, i);
                return m_data.set_data(pair.first, pair.second ? consts::conj(elem) : elem);
            }
        }

        T get(size_t a, size_t i) const override {
            // any compiler should statically execute this conditional
            if (!consts::is_complex<T>())
                return m_data.get_data(m_indexer.index_only(a, i));
            else {
                const auto pair = m_indexer.index_and_conj(a, i);
                const auto element = m_data.get_data(pair.first);
                return pair.second ? consts::conj(element) : element;
            }
        }
    };

    template<typename T> using SymNone = IndexedArray<IndexerSymNone, T>;
    template<typename T> using SymH = IndexedArray<IndexerSymH, T>;

}

#endif //M7_INTEGRALARRAY1E_H
