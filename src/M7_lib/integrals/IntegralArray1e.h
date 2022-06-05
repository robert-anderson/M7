//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY1E_H
#define M7_INTEGRALARRAY1E_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/defs.h>
#include "IntegralStorage.h"

namespace integrals_1e {

    struct Indexer : IntegralIndexer {
        Indexer(size_t norb, size_t size): IntegralIndexer(norb, size){}
        virtual void foreach(const std::function<void(size_t, size_t)>& fn) const = 0;
    };

    struct IndexerSymNone : Indexer {
        IndexerSymNone(size_t norb);

        size_t index_only(size_t a, size_t i) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t i) const;

        void foreach(const std::function<void(size_t, size_t)> &fn) const override;
    };

    using namespace integer_utils;

    struct IndexerSymH : Indexer {
        IndexerSymH(size_t norb);

        size_t index_only(size_t a, size_t i) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t i) const;

        void foreach(const std::function<void(size_t, size_t)> &fn) const override;
    };

    template<typename T>
    struct Array {
        const size_t m_norb;
        Array(size_t norb): m_norb(norb){}
        virtual ~Array() = default;

        virtual bool set(size_t a, size_t i, T elem) = 0;

        virtual T get(size_t a, size_t i) const = 0;

        virtual void transfer(const Array<T>* higher_sym) = 0;

    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<Indexer, indexer_t>::value, "invalid template class");
        using Array<T>::set;
        using Array<T>::get;
        indexer_t m_indexer;
        PrivateIntegralStorage<T> m_data;

        IndexedArray(size_t norb) : Array<T>(norb), m_indexer(norb),
            m_data(static_cast<const Indexer&>(m_indexer).m_size) {}

        void transfer(const Array<T>* higher_sym) override {
            auto fn = [&](size_t a, size_t i){
                bool success = set(a, i, higher_sym->get(a, i));
                REQUIRE_TRUE(success, "error in transferring 2e integral array");
            };
            static_cast<const Indexer&>(m_indexer).foreach(fn);
        }

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

        std::string name(Sym sym);
    }

    /**
     * @tparam T
     *  matrix element type
     * @param norb
     *  number of orbitals to pass to ctor
     * @param sym
     *  the symmetry to use in the new allocation (decremented at output)
     * @return
     *  new instance of the type corresponding to the requested symmetry
     */
    template<typename T>
    std::unique_ptr<Array<T>> make(size_t norb, syms::Sym sym){
        typedef std::unique_ptr<Array<T>> ptr_t;
        using namespace syms;
        switch (sym) {
            case(Null):
                return {};
            case(None):
                return ptr_t(new SymNone<T>(norb));
            case(H):
                return ptr_t(new SymH<T>(norb));
            default:
                return {};
        }
    }

    /**
     * when attempting to fill the integral arrays, the highest symmetry is assumed at first. when a counter example to
     * that symmetry assumption is reached, the parse is restarted with the next lower symmetry until the entire
     * integral array can be stored without contradiction
     * @tparam T
     *  matrix element type
     * @param ptr
     *  smart pointer to the currently-allocated array
     * @param sym
     *  the symmetry to use in the new allocation (decremented at output)
     */
    template<typename T>
    void next_sym_attempt(std::unique_ptr<Array<T>>& ptr, syms::Sym& sym){
        if (!ptr) return;
        std::unique_ptr<Array<T>> new_ptr = make<T>(ptr->m_norb, sym);
        if (!new_ptr) {
            ptr = nullptr;
            return;
        }
        new_ptr->transfer(ptr.get());
        ptr = std::move(new_ptr);
        sym = syms::Sym(sym-1);
    }
}

#endif //M7_INTEGRALARRAY1E_H
