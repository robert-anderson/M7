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
 * KHDRDR = 1
 * KHD = RDR
 *
 * HD = K <ji|ba>
 * KHD = RDR
 *
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
 *  [AB, C] = [A, B]C + B[A, C]
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

    struct IndexerSymNone : IntegralIndexer {
        const size_t m_norb2, m_norb3;

        IndexerSymNone(size_t norb) : IntegralIndexer(norb, pow(norb, 4ul)),
                                      m_norb2(norb * norb), m_norb3(norb * m_norb2) {}

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            return a * m_norb3 + b * m_norb2 + i * m_norb + j;
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            return {index_only(a, b, i, j), false};
        }
    };

    struct IndexerSymH : IntegralIndexer {
        IndexerSymH(size_t norb) : IntegralIndexer(norb, npair(norb*norb)) {}

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            const auto ab = a*m_norb+b;
            const auto ij = i*m_norb+j;
            return (ab>=ij) ? trigmap(ab, ij) : trigmap(ij, ab);
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            const auto ab = a*m_norb+b;
            const auto ij = i*m_norb+j;
            if (ab>=ij) return {trigmap(ab, ij), true};
            else return {trigmap(ij, ab), false};
        }
    };

    struct IndexerSymD : IntegralIndexer {
        IndexerSymH m_sym_h;
        IndexerSymD(size_t norb) : IntegralIndexer(norb, npair(norb*norb)), m_sym_h(norb){}

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            return m_sym_h.index_only(a, i, b, j);
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            return {index_only(a, b, i, j), false};
        }
    };

    struct IndexerSymDHR : IntegralIndexer {
        IndexerSymDHR(size_t norb) : IntegralIndexer(norb, npair(npair(norb))){}

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            auto ab = (a>=b) ? trigmap(a, b) : trigmap(b, a);
            auto ij = (i>=j) ? trigmap(i, j) : trigmap(j, i);
            return (ab>=ij) ? trigmap(ab, ij) : trigmap(ij, ab);
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            // always have real orbitals so never any conjugation
            return {index_only(a, b, i, j), false};
        }
    };

    /**
     * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays.
     * in the HI case, we have one eightfold array corresponding to elements equivalent to:
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
    struct IndexerSymDH : IntegralIndexer {
        const size_t m_hir_size;
        IndexerSymDH(size_t norb) :  IntegralIndexer(norb, 2*npair(npair(norb))), m_hir_size(m_size/2){}

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            // let y = <ab|ij>, correct order is x
            if (b >= j){
                const auto bj = trigmap(b, j);
                if (a >= i) {
                    const auto ai= trigmap(a, i);
                    if (bj >= ai) {
                        // x = y
                        return {trigmap(bj, ai), false};
                    }
                    else {
                        // x = Dy i.e. y = Dx
                        return {trigmap(ai, bj), false};
                    }
                }
                else {
                    const auto ai= trigmap(i, a);
                    if (bj >= ai) {
                        // x = KHRy i.e. y = RHKx = HRKx
                        return {trigmap(bj, ai)+m_hir_size, true};
                    }
                    else {
                        // x = DKHRy i.e. y = RHDKx = DRx
                        return {trigmap(ai, bj)+m_hir_size, false};
                    }
                }
            }
            else {
                const auto bj = trigmap(j, b);
                if (a >= i){
                    const auto ai= trigmap(a, i);
                    if (bj>=ai) {
                        // x = Ry i.e. y = Rx
                        return {trigmap(bj, ai)+m_hir_size, false};
                    }
                    else {
                        // x = DRy i.e. y = RDx = HDRKx
                        return {trigmap(ai, bj)+m_hir_size, true};
                    }
                }
                else {
                    const auto ai= trigmap(i, a);
                    if (bj>=ai) {
                        // x = KHRRy i.e. y = HKx
                        return {trigmap(bj, ai), true};
                    }
                    else {
                        // x = DKHRRy i.e. y = HDKx
                        return {trigmap(ai, bj), true};
                    }
                }
            }
            return {};
        }

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            return index_and_conj(a, b, i, j).first;
        }

    };

    /**
     * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays.
     * in the DR case, we have one eightfold array corresponding to elements equivalent to:
     *      (x)             (Rx)            (Dx)          (DRx)
     *  <ab|ij>  =  <aj|ib>  =  <ba|ji>  =  <ja|bi>
     *
     *  and another corresponding to elements equivalent to:
     *      (Hx)            (RHx)          (DHx)         (DRHx)
     *  <ij|ab>  =  <ib|aj>  =  <ji|ba>  =  <bi|ja>
     *
     *  determining the index (data offset) and whether to apply complex conjugation is simply a case of keeping track
     *  of which of the permuting operations (H, I, R) are required to bring the given indices to the canonical ordering
     *  i.e. a>=b, i>=j, ab>=ij
     */
    struct IndexerSymDR : IntegralIndexer {
        const size_t m_hir_size;

        IndexerSymDR(size_t norb) : IntegralIndexer(norb, 2 * npair(npair(norb))), m_hir_size(m_size / 2) {}

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            return {};
        }

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            return index_and_conj(a, b, i, j).first;
        }

    };

#if 0
    struct IndexerSymIR : IntegbralIndexer {
        IndexerSymHI m_sym_hi;
        IndexerSymIR(size_t norb) :  IntegralIndexer(norb, 2*npair(npair(norb))), m_sym_hi(norb){}

        size_t index_only(size_t a, size_t b, size_t i, size_t j) const {
            return m_sym_hi.index_only(a, i, b, j);
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t b, size_t i, size_t j) const {
            return {};
        }
    };
#endif



#if 0



    using namespace integer_utils;
    struct Integral1eIndexerH : Integral1eIndexer {
        Integral1eIndexerH(size_t norb): Integral1eIndexer(norb, npair(norb)){}

        size_t get_case(size_t a, size_t i) const {
            return a >= i;
        }

        size_t index_only(size_t a, size_t i) const {
            switch (get_case(a, i)) {
                case 0: return trigmap(i, a);
                case 1: return trigmap(a, i);
            }
            return ~0ul;
        }

        std::pair<size_t, bool> index_and_conj(size_t a, size_t i) const {
            switch (get_case(a, i)) {
                case 0: return {trigmap(i, a), true};
                case 1: return {trigmap(a, i), false};
            }
            return {~0ul, false};
        }
    };
#endif


    template<typename T>
    struct Array {
        virtual bool set(size_t a, size_t b, size_t i, size_t j, T elem) = 0;

        virtual T get(size_t a, size_t b, size_t i, size_t j) const = 0;
    };

    template<typename indexer_t, typename T>
    struct IndexedArray : Array<T> {
        static_assert(std::is_base_of<IntegralIndexer, indexer_t>::value, "invalid template class");
        using Array<T>::set;
        using Array<T>::get;
        indexer_t m_indexer;
        SharedIntegralStorage<T> m_data;

        IndexedArray(size_t norb) : m_indexer(norb), m_data(static_cast<const IntegralIndexer &>(m_indexer).m_size) {}

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
    };

    template<typename T> using SymNone = IndexedArray<IndexerSymNone, T>;
    template<typename T> using SymH = IndexedArray<IndexerSymH, T>;
    template<typename T> using SymI = IndexedArray<IndexerSymD, T>;
    template<typename T> using SymDH = IndexedArray<IndexerSymDH, T>;
    template<typename T> using SymDR = IndexedArray<IndexerSymDR, T>;
    template<typename T> using SymHIR = IndexedArray<IndexerSymDHR, T>;
}

#if 0

const size_t m_n8fold;

private:





/**
 * store two-electron integrals with the physicists' notation (ab|ij) where the integrand is:
 *      a*(x1)b*(x2) 1/r12 i(x1)j(x2)
 * such quantities admit varying levels of permutational symmetry whereby an integral indexed by <ab|ij> is equal (up to
 * a complex conjugation) to any of the following
 *
 *  <ba|ji>     (fermion indistinguishability, or integration variable interchange "i")
 *  <ij|ab>     (hermiticity "h")
 *  <ab|ji>     (real orbitals "R")
 *
 * The assumption of:
 *  - of "i" only is called "2-fold" symmetry,
 *  - of "ih" is a "4-fold" symmetry for complex orbitals
 *  - of "ir" is another "4-fold" symmetry for non-hermitian hamiltonians
 *  - and of "ihr" is called "8-fold symmetry"
 *
 *  "i" computes combined ai and jb indices, and only stores the elements with ai>=jb
 *
 *
 *  8 fold only stores the elements with i>=j and k>=l and ij>=kl
 *  4 fold is a slightly more complicated case, which is important in situations where the integrals correspond to
 *      complex-valued orbitals. the full explanation is given in the subclass docs
 */
template<typename T>
struct IntegralArray2e : IntegralStorage<T> {
protected:
    using typename IntegralStorage<T>::real_t;
    using typename IntegralStorage<T>::cmplx_t;
    typedef SharedArray<real_t> real_array_t;
    typedef SharedArray<cmplx_t> cmplx_array_t;
    SharedArray<T> m_data;
    IntegralArray2e(size_t norb, size_t size): IntegralStorage<T>(norb), m_data(size){}

private:
    
    virtual void ref_get(size_t a, size_t b, size_t i, size_t j, const real_array_t &data, real_t& elem) const = 0;
    virtual void ref_get(size_t a, size_t b, size_t i, size_t j, const cmplx_array_t &data, cmplx_t& elem) const = 0;
    virtual void ref_set(size_t a, size_t b, size_t i, size_t j, real_array_t &data, const real_t& elem) const = 0;
    virtual void ref_set(size_t a, size_t b, size_t i, size_t j, cmplx_array_t &data, const cmplx_t& elem) const = 0;

public:
    T get(size_t a, size_t b, size_t i, size_t j) const {
        T elem;
        ref_get(a, b, i, j, m_data, elem);
        return elem;
    }

    void set(size_t a, size_t b, size_t i, size_t j, T elem) {
        ref_set(a, b, i, j, m_data, elem);
    }
};


template<typename T>
struct IntegralArray2e_i : IntegralArray2e<T> {
    using IntegralArrayBase::trig;
    using IntegralArrayBase::m_norb;
    using typename IntegralStorage<T>::real_t;
    using typename IntegralStorage<T>::cmplx_t;
    using typename IntegralArray2e<T>::real_array_t;
    using typename IntegralArray2e<T>::cmplx_array_t;
    IntegralArray2e_i(size_t norb) : IntegralArray2e<T>(norb, trig(norb*norb, 0)){}

private:
    size_t flatten(size_t a, size_t b, size_t i, size_t j) const {
        const auto ai = a*m_norb+i;
        const auto bj = b*m_norb+j;
        return ai>=bj ? trig(ai, bj) : trig(bj, ai);
    }

    void ref_get(size_t a, size_t b, size_t i, size_t j,
                 const real_array_t &data, real_t &elem) const override {
        elem = data[flatten(a, b, i, j)];
    }

    void ref_get(size_t a, size_t b, size_t i, size_t j, const cmplx_array_t &data, cmplx_t &elem) const override {
        elem = data[flatten(a, b, i, j)];
    }

    void ref_set(size_t a, size_t b, size_t i, size_t j, real_array_t &data, const real_t &elem) const override {
        data.set(flatten(a, b, i, j), elem);
    }

    void ref_set(size_t a, size_t b, size_t i, size_t j, cmplx_array_t &data, const cmplx_t &elem) const override {
        data.set(flatten(a, b, i, j), elem);
    }
};

template<typename T>
struct IntegralArray2e_ihr : IntegralArray2e<T> {
    using IntegralArrayBase::trig;
    using typename IntegralStorage<T>::real_t;
    using typename IntegralStorage<T>::cmplx_t;
    using typename IntegralArray2e<T>::real_array_t;
    using typename IntegralArray2e<T>::cmplx_array_t;
    IntegralArray2e_ihr(size_t norb): IntegralArray2e<T>(norb, trig(trig(norb, 0), trig(norb, 0))){}

private:

    static size_t flatten(size_t a, size_t b, size_t i, size_t j) {
        if (a >= b){
            const auto ab = trig(a, b);
            if (i >= j){
                const auto ij = trig(i, j);
                //                  <ab|ij>        <ij|ab>
                return (ab >= ij) ? trig(ab, ij) : trig(ij, ab);
            }
            else {
                auto ij = trig(j, i);
                //                  <ab|ji>        <ji|ab>
                return (ab >= ij) ? trig(ab, ij) : trig(ij, ab);
            }
        }
        else {
            auto ab = trig(b, a);
            if (i >= j){
                auto ij = trig(i, j);
                //                  <ba|ij>        <ij|ba>
                return (ab >= ij) ? trig(ab, ij) : trig(ij, ab);
            }
            else {
                auto ij = trig(j, i);
                //                  <ba|ji>        <ji|ba>
                return (ab >= ij) ? trig(ab, ij) : trig(ij, ab);
            }
        }
        return ~0ul;
    }


    void ref_get(size_t a, size_t b, size_t i, size_t j, const real_array_t &data, real_t &elem) const override {
        elem = data[flatten(a, b, i, j)];
    }

    /**
     * it's impossible for real orbitals to yield integrals with non-zero imaginary part, so in this case it must be
     * assumed that we are (wastefully) using a complex container for real-valued integrals
     */
    void ref_get(size_t a, size_t b, size_t i, size_t j, const cmplx_array_t &data, cmplx_t &elem) const override {
        elem = data[flatten(a, b, i, j)];
    }

    void ref_set(size_t a, size_t b, size_t i, size_t j, real_array_t &data, const real_t &elem) const override {
        data.set(flatten(a, b, i, j), elem);
    }

    void ref_set(size_t a, size_t b, size_t i, size_t j, cmplx_array_t &data, const cmplx_t &elem) const override {
        data.set(flatten(a, b, i, j), elem);
    }
};

/**
 * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays,
 *                                                    "I"     "H"          "IH"
 * one for the indices equivalent to (ij|kl), i.e. (kl|ij), (ji|lk), and (lk|ji)
 * and another for those equivalent to (ij|lk) i.e. (lk|ij), (ji|kl), and (kl|ji)
 *
 *                                                    "I"     "R"          "RI"
 * one for the indices equivalent to (ij|kl), i.e. (kl|ij), (ij|lk), and (kl|ji)
 * and another for those equivalent to (ji|lk) i.e. (lk|ji), (ji|kl), and (lk|ij)
 *
 *
 */

template<typename T>
struct IntegralArray2e_ih : IntegralArray2e<T> {
    using IntegralArrayBase::trig;
    using typename IntegralStorage<T>::real_t;
    using typename IntegralStorage<T>::cmplx_t;
    using typename IntegralArray2e<T>::real_array_t;
    using typename IntegralArray2e<T>::cmplx_array_t;
    using IntegralArray2e<T>::m_data;
    const size_t m_n8fold;
    IntegralArray2e_ih(size_t norb):
        IntegralArray2e<T>(norb, 2*trig(trig(norb, 0), trig(norb, 0))),
        m_n8fold(m_data.size()/2){}

private:

    size_t flatten(size_t i, size_t j, size_t k, size_t l) const {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                return (ij>=kl) ? trig(ij, kl) : trig(kl, ij);
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)                        (lk|ij)
                return ((ij>=kl) ? trig(ij, kl) : trig(kl, ij)) + m_n8fold;
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                     (ji|kl)         (kl|ji)
                return ((ij>=kl) ? trig(ij, kl) : trig(kl, ij)) + m_n8fold;
            }
            else {
                auto kl = trig(l, k);
                //                    (ji|lk)        (lk|ji)
                return (ij>=kl) ? trig(ij, kl): trig(kl, ij);
            }
        }
        return ~0ul;
    }

    void ref_get(size_t i, size_t j, size_t k, size_t l, const real_array_t &data, real_t &elem) const override {
        elem = data[flatten(i, j, k, l)];
    }

    void ref_get(size_t i, size_t j, size_t k, size_t l, const cmplx_array_t &data, cmplx_t &elem) const override {
        elem = data[flatten(i, j, k, l)];
    }

    void ref_set(size_t i, size_t j, size_t k, size_t l, real_array_t &data, const real_t &elem) const override {
        data.set(flatten(i, j, k, l), elem);
    }

    void ref_set(size_t i, size_t j, size_t k, size_t l, cmplx_array_t &data, const cmplx_t &elem) const override {
        data.set(flatten(i, j, k, l), elem);
    }
};


#endif //M7_INTEGRALARRAY2E_H
#endif //M7_INTEGRALARRAY2E_H