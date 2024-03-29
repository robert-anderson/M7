//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY2E_H
#define M7_INTEGRALARRAY2E_H

#include <M7_lib/defs.h>
#include "IntegralStorage.h"
#include <M7_lib/util/Integer.h>

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
    using namespace integer;

    typedef std::function<void(uint_t, uint_t, uint_t, uint_t)> foreach_fn_t;

    namespace syms {
        enum Sym {
            Null, None, H, D, DR, DH, DHR
        };

        str_t desc(Sym sym);

        Sym from_symbol(str_t name);

        strv_t equivalences(Sym sym);
    }

    struct Indexer : IntegralIndexer {
        const syms::Sym m_sym;
        Indexer(uint_t norb, uint_t size, syms::Sym sym): IntegralIndexer(norb, size), m_sym(sym){}

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
        const uint_t m_norb2, m_norb3;

        IndexerSymNone(uint_t norb);

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;

    };

    /**
     * indexing scheme which only assumes hermiticity. this is another unphysical scheme for the same reason as above
     */
    struct IndexerSymH : Indexer {
        IndexerSymH(uint_t norb);

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;

    };

    /**
     * lowest-symmetry conceivably useful indexing scheme. for use when orbitals are complex, and H is non-hermitian
     */
    struct IndexerSymD : Indexer {
        IndexerSymH m_sym_h;
        IndexerSymD(uint_t norb);

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;
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
        const uint_t m_hir_size;
        IndexerSymDH(uint_t norb);

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

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
        const uint_t m_hir_size;
        IndexerSymDR(uint_t norb);

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;
    };

    /**
     * indexer with full "8-fold" symmetry. for use with most non-relativistic ab-initio Hamiltonians:
     *  - real orbitals
     *  - Hermitian H
     */
    struct IndexerSymDHR : Indexer {
        IndexerSymDHR(uint_t norb);

        uint_t index_only(uint_t a, uint_t b, uint_t i, uint_t j) const;

        std::pair<uint_t, bool> index_and_conj(uint_t a, uint_t b, uint_t i, uint_t j) const;
    };


    template<typename T>
    struct Array {
        const uint_t m_norb;
        Array(uint_t norb): m_norb(norb){}
        virtual ~Array() = default;

        virtual bool set_(uint_t a, uint_t b, uint_t i, uint_t j, T elem) = 0;

        virtual T get(uint_t a, uint_t b, uint_t i, uint_t j) const = 0;

        virtual syms::Sym sym() const = 0;
    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<IntegralIndexer, indexer_t>::value, "invalid template class");
        using Array<T>::set_;
        using Array<T>::get;
        indexer_t m_indexer;
        SharedIntegralStorage<T> m_data;

        IndexedArray(uint_t norb) : Array<T>(norb), m_indexer(norb), m_data(m_indexer.m_size) {}

        bool set_(uint_t a, uint_t b, uint_t i, uint_t j, T elem) override {
            // any compiler should statically execute this conditional
            if (!dtype::is_complex<T>())
                return m_data.set_data_(m_indexer.index_only(a, b, i, j), elem);
            else {
                const auto pair = m_indexer.index_and_conj(a, b, i, j);
                return m_data.set_data_(pair.first, pair.second ? arith::conj(elem) : elem);
            }
        }

        T get(uint_t a, uint_t b, uint_t i, uint_t j) const override {
            // any compiler should statically execute this conditional
            if (!dtype::is_complex<T>())
                return m_data.get_data(m_indexer.index_only(a, b, i, j));
            else {
                const auto pair = m_indexer.index_and_conj(a, b, i, j);
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
    template<typename T> using SymD = IndexedArray<IndexerSymD, T>;
    template<typename T> using SymDH = IndexedArray<IndexerSymDH, T>;
    template<typename T> using SymDR = IndexedArray<IndexerSymDR, T>;
    template<typename T> using SymDHR = IndexedArray<IndexerSymDHR, T>;

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
            case D:
                return ptr_t(new SymD<T>(norb));
            case DH:
                return ptr_t(new SymDH<T>(norb));
            case DR:
                return ptr_t(new SymDR<T>(norb));
            case DHR:
                return ptr_t(new SymDHR<T>(norb));
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
        REQUIRE_FALSE(ptr->sym()==syms::None, "Integrals are internally inconsistent");
        if (ptr->sym()==syms::D)
            logging::warn("Integration variable interchange should be a symmetry of the 2-electron integrals!");
        std::unique_ptr<Array<T>> new_ptr = make<T>(ptr->m_norb, syms::Sym(ptr->sym()-1));
        if (!new_ptr) {
            ptr = nullptr;
            return;
        }
        ptr = std::move(new_ptr);
    }
}

#endif //M7_INTEGRALARRAY2E_H