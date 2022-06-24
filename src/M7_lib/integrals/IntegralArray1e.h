//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY1E_H
#define M7_INTEGRALARRAY1E_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/defs.h>
#include "IntegralStorage.h"
#include <M7_lib/util/Integer.h>

namespace integrals_1e {
    using namespace integer;

    typedef std::function<void(uint_t, uint_t)> foreach_fn_t;

    namespace syms {
        enum Sym {
            Null, None, H
        };

        std::string name(Sym sym);
    }


    struct Indexer : IntegralIndexer {
        const syms::Sym m_sym;
        Indexer(uint_t norb, uint_t size, syms::Sym sym): IntegralIndexer(norb, size), m_sym(sym){}
        virtual void foreach(const foreach_fn_t& fn) const = 0;
    };

    struct IndexerSymNone : Indexer {
        IndexerSymNone(uint_t norb);

        uint_t index_only(uint_t a, uint_t i) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t i) const;

        void foreach(const foreach_fn_t &fn) const override;
    };

    using namespace integer;

    struct IndexerSymH : Indexer {
        IndexerSymH(uint_t norb);

        uint_t index_only(uint_t a, uint_t i) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t i) const;

        void foreach(const foreach_fn_t &fn) const override;
    };

    template<typename T>
    struct Array {
        const uint_t m_norb;
        Array(uint_t norb): m_norb(norb){}
        virtual ~Array() = default;

        virtual bool set(uint_t a, uint_t i, T elem) = 0;

        virtual T get(uint_t a, uint_t i) const = 0;

        virtual syms::Sym sym() const = 0;
    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<Indexer, indexer_t>::value, "invalid template class");
        using Array<T>::set;
        using Array<T>::get;
        indexer_t m_indexer;
        PrivateIntegralStorage<T> m_data;

        IndexedArray(uint_t norb) : Array<T>(norb), m_indexer(norb),
            m_data(static_cast<const Indexer&>(m_indexer).m_size) {}

        bool set(uint_t a, uint_t i, T elem) override {
            // any compiler should statically execute this conditional
            if (!dtype::is_complex<T>())
                return m_data.set_data(m_indexer.index_only(a, i), elem);
            else {
                const auto pair = m_indexer.index_and_conj(a, i);
                return m_data.set_data(pair.first, pair.second ? arith::conj(elem) : elem);
            }
        }

        T get(uint_t a, uint_t i) const override {
            // any compiler should statically execute this conditional
            if (!dtype::is_complex<T>())
                return m_data.get_data(m_indexer.index_only(a, i));
            else {
                const auto pair = m_indexer.index_and_conj(a, i);
                const auto element = m_data.get_data(pair.first);
                return pair.second ? arith::conj(element) : element;
            }
        }

        syms::Sym sym() const override {
            return static_cast<const Indexer&>(m_indexer).m_sym;
        }

    };

    template<typename T> using SymNone = IndexedArray<IndexerSymNone, T>;
    template<typename T> using SymH = IndexedArray<IndexerSymH, T>;

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
    std::unique_ptr<Array<T>> make(uint_t norb, syms::Sym sym){
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
     */
    template<typename T>
    void next_sym_attempt(std::unique_ptr<Array<T>>& ptr){
        if (!ptr) return;
        std::unique_ptr<Array<T>> new_ptr = make<T>(ptr->m_norb, syms::Sym(ptr->sym()-1));
        if (!new_ptr) {
            ptr = nullptr;
            return;
        }
        ptr = std::move(new_ptr);
    }
}

#endif //M7_INTEGRALARRAY1E_H
