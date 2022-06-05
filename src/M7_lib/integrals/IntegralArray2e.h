//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY2E_H
#define M7_INTEGRALARRAY2E_H

#include <M7_lib/defs.h>
#include "IntegralStorage.h"

/**
 * all two-electron integrals are denoted in physicists' ordering: <ab|ij>
 * there are three permutational operators associated with symmetry in these objects:
 *  - D: interchange of dummy variables of integration I <ab|ij> = <ba|ji>
 *  - H: hermitian transposition H <ab|ij> = K <ij|ab> where K is the complex conjugation operator
 *  - R: conjugation of one-electron functions of one dummy variable R <ab|ij> = <aj|ib>. R is chosen to act
 *       on the second pair of functions, so KHR acts on the first
 * all operators are their own inverses
 *
 * as shown below, all commutators are zero with the exception of [R, D].
 *
 * HR <ab|ij> = H <aj|ib> = K <ib|aj>
 * DR <ab|ij> = D <aj|ib> = <ja|bi>
 * RD <ab|ij> = R <ba|ji> = <bi|ja>
 * DRD <ab|ij> = <ib|aj>
 * KHR <ab|ij> = <ib|aj>
 *
 * i.e. [R, D] != 0
 * RDR <ab|ij> = <ji|ba>
 * DRDR <ab|ij> = <ij|ab>
 * HDRDR <ab|ij> = K <ab|ij>
 *
 * i.e. KHDRDR = 1
 *
 * K <ab|ij> = H <ij|ab>
 * HK <ab|ij> = HH <ij|ab> = <ij|ab>
 * KH <ab|ij> = KK<ij|ab> = <ij|ab>
 * so [H, K] = 0
 *
 * DKH <ab|ij> = D <ij|ab> = <ji|ba>
 * KHD <ab|ij> = KH <ba|ji> = <ji|ba>
 *
 * i.e. [KH, D] = 0
 *  in general [AB, C] = [A, B]C + B[A, C]
 *  [KH, D] = [K, H]D + H[K, D]
 *  0 = 0 + H[K, D]
 * i.e. [K, D] = 0
 *
 * KHR <ab|ij> = KH <aj|ib> = <ib|aj>
 * RKH <ab|ij> = R <ij|ab> = <ib|aj>
 *
 * i.e. [KH, R] = 0
 *  [KH, R] = [K, H]R + H[K, R]
 * i.e. [K, R] = [R, H] = 0
 *
 */
namespace integrals_2e {
    using namespace integer_utils;

    typedef std::function<void(size_t, size_t, size_t, size_t)> foreach_fn_t;

    struct Indexer : IntegralIndexer {
        Indexer(size_t norb, size_t size): IntegralIndexer(norb, size){}

        /**
         * currently this is just a general (assuming no symmetry) iterator, and could be overridden in derived classes,
         * but since this is not currently used in any performance-critical sector of the code, these symmetry-aware
         * overloads are not implemented
         * @param fn
         *  function to execute on each quadruplet of integral indices
         */
        virtual void foreach(const foreach_fn_t& fn) const;
    };

    /**
     * basic (practically unused) indexing scheme which assumes no permutational symmetries.
     * this is actually an unphysical assumption, since any two-electron integrals will have D symmetry to exploit
     */
    struct IndexerSymNone : Indexer {
        const size_t m_norb2, m_norb3;

        IndexerSymNone(size_t norb);

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;

    };

    /**
     * indexing scheme which only assumes hermiticity. this is another unphysical scheme for the same reason as above
     */
    struct IndexerSymH : Indexer {
        IndexerSymH(size_t norb);

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;

    };

    /**
     * lowest-symmetry conceivably useful indexing scheme. for use when orbitals are complex, and H is non-hermitian
     */
    struct IndexerSymD : Indexer {
        IndexerSymH m_sym_h;
        IndexerSymD(size_t norb);

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;
    };

    /**
     * Indexing scheme which does not assume real orbitals, but assumes hermiticity.
     * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays.
     * in the DH case, we have one eightfold array corresponding to elements equivalent to:
     *    (x)         (Dx)       (HKx)      (HDKx)
     *  <ab|ij>  =  <ba|ji>  =  <ij|ab>  =  <ji|ba>
     *
     *  and another corresponding to elements equivalent to:
     *    (Rx)       (DRx)      (HRKx)      (HDRKx)
     *  <aj|ib>  =  <ja|bi>  =  <ib|aj>  =  <bi|ja>
     *
     *  we have to fuse index pairs which can be reversed under the application of R and DR: (bj, ai)
     *
     *  determining the index (data offset) and whether to apply complex conjugation is simply a case of keeping track
     *  of which of the permuting operations (H, I, R) are required to bring the given indices to the canonical ordering
     *  i.e. b>=j, a>=i, bj>=ai
     */
    struct IndexerSymDH : Indexer {
        const size_t m_hir_size;
        IndexerSymDH(size_t norb);

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

    };

    /**
     * Indexing scheme which does not assume hermiticity, but assumes real orbitals.
     * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays.
     * in the DR case, we have one eightfold array corresponding to elements equivalent to:
     *    (x)         (Rx)       (Dx)        (DRx)
     *  <ab|ij>  =  <aj|ib>  =  <ba|ji>  =  <ja|bi>
     *
     *  and another corresponding to elements equivalent to:
     *   (HKx)      (HRKx)      (HDKx)      (HDRKx)
     *  <ij|ab>  =  <ib|aj>  =  <ji|ba>  =  <bi|ja>
     *
     *  we have to fuse index pairs which can be reversed under the application of R and DR: (bj, ai)
     *
     *  determining the index (data offset) and whether to apply complex conjugation is simply a case of keeping track
     *  of which of the permuting operations (H, I, R) are required to bring the given indices to the canonical ordering
     *  i.e. b>=j, a>=i, bj>=ai
     */
    struct IndexerSymDR : Indexer {
        const size_t m_hir_size;
        IndexerSymDR(size_t norb);

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;
    };

    /**
     * indexer with full "8-fold" symmetry. for use with most non-relativistic ab-initio Hamiltonians:
     *  - real orbitals
     *  - Hermitian H
     */
    struct IndexerSymDHR : Indexer {
        IndexerSymDHR(size_t norb);

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const;

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const;
    };


    template<typename T>
    struct Array {
        const size_t m_norb;
        Array(size_t norb): m_norb(norb){}
        virtual ~Array() = default;

        virtual bool set(size_t a, size_t b, size_t i, size_t j, T elem) = 0;

        virtual T get(size_t a, size_t b, size_t i, size_t j) const = 0;

        virtual void transfer(const Array<T>* higher_sym) = 0;
    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<IntegralIndexer, indexer_t>::value, "invalid template class");
        using Array<T>::set;
        using Array<T>::get;
        indexer_t m_indexer;
        SharedIntegralStorage<T> m_data;

        IndexedArray(size_t norb) :
            Array<T>(norb), m_indexer(norb), m_data(static_cast<const IntegralIndexer &>(m_indexer).m_size) {}

        bool set(size_t a, size_t b, size_t i, size_t j, T elem) override {
            // any compiler should statically execute this conditional
            if (!consts::is_complex<T>())
                return m_data.set_data(m_indexer.index_only(a, b, i, j), elem);
            else {
                const auto pair = m_indexer.index_and_conj(a, b, i, j);
                return m_data.set_data(pair.first, pair.second ? consts::conj(elem) : elem);
            }
        }

        T get(size_t a, size_t b, size_t i, size_t j) const override {
            // any compiler should statically execute this conditional
            if (!consts::is_complex<T>())
                return m_data.get_data(m_indexer.index_only(a, b, i, j));
            else {
                const auto pair = m_indexer.index_and_conj(a, b, i, j);
                const auto element = m_data.get_data(pair.first);
                return pair.second ? consts::conj(element) : element;
            }
        }

        void transfer(const Array<T> *higher_sym) override {
            auto fn = [&](size_t a, size_t b, size_t i, size_t j) {
                bool success = set(a, b, i, j, higher_sym->get(a, b, i, j));
                REQUIRE_TRUE(success, "error in transferring integral array");
            };
            static_cast<const Indexer&>(m_indexer).foreach(fn);
        }
    };

    template<typename T> using SymNone = IndexedArray<IndexerSymNone, T>;
    template<typename T> using SymH = IndexedArray<IndexerSymH, T>;
    template<typename T> using SymD = IndexedArray<IndexerSymD, T>;
    template<typename T> using SymDH = IndexedArray<IndexerSymDH, T>;
    template<typename T> using SymDR = IndexedArray<IndexerSymDR, T>;
    template<typename T> using SymDHR = IndexedArray<IndexerSymDHR, T>;

    namespace syms {
        enum Sym {
            Null, None, H, D, DH, DR, DHR
        };

        std::string name(Sym sym);

        std::vector<std::string> equivalences(Sym sym);
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
        typedef std::unique_ptr<Array<T>> ptr_t;
        ptr_t new_ptr = nullptr;
        using namespace syms;
        const auto norb = ptr->m_norb;
        switch (sym) {
            case(Null):
                ptr = nullptr;
                return;
            case(None):
                new_ptr = ptr_t(new SymNone<defs::ham_t>(norb));
                break;
            case(H):
                new_ptr = ptr_t(new SymH<defs::ham_t>(norb));
                break;
            case D:
                new_ptr = ptr_t(new SymD<defs::ham_t>(norb));
                break;
            case DH:
                new_ptr = ptr_t(new SymDH<defs::ham_t>(norb));
                break;
            case DR:
                new_ptr = ptr_t(new SymDR<defs::ham_t>(norb));
                break;
            case DHR:
                new_ptr = ptr_t(new SymDHR<defs::ham_t>(norb));
                break;
        }
        new_ptr->transfer(ptr.get());
        ptr = std::move(new_ptr);
        sym = Sym(sym-1);
    }

}

#endif //M7_INTEGRALARRAY2E_H