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
        Array() = default;
        virtual ~Array() = default;

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

    namespace syms {
        enum Sym {
            Null, None, H
        };
    }

    /**
     * when attempting to fill the integral arrays, the highest symmetry is assumed at first. when a counter example to
     * that symmetry assumption is reached, the parse is restarted with the next lower symmetry until the entire
     * integral array can be stored without contradiction
     * @tparam T
     *  matrix element type
     * @param norb
     *  number of one-electron (spinorbitals if spin resolved basis) functions in each dimension of the integral arrays
     * @param ptr
     *  smart pointer to the currently-allocated array
     * @param sym
     *  the symmetry to use in the new allocation (decremented at output)
     */
    template<typename T>
    void next_sym_attempt(size_t norb, std::unique_ptr<Array<T>>& ptr, syms::Sym sym){
        typedef std::unique_ptr<Array<T>> ptr_t;
        using namespace syms;
        switch (sym) {
            case(Null):
                ptr = nullptr;
                return;
            case(None):
                ptr = ptr_t(new SymNone<T>(norb));
            case(H):
                ptr = ptr_t(new SymH<T>(norb));
        }
        sym = Sym(sym-1);
    }
}

#endif //M7_INTEGRALARRAY1E_H
