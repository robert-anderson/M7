//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY2E_H
#define M7_INTEGRALARRAY2E_H

#include <M7_lib/defs.h>
#include <M7_lib/parallel/SharedArray.h>

#include "IntegralArray.h"

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
struct IntegralArray2e : IntegralArray<T> {
protected:
    using typename IntegralArray<T>::real_t;
    using typename IntegralArray<T>::cmplx_t;
    typedef SharedArray<real_t> real_array_t;
    typedef SharedArray<cmplx_t> cmplx_array_t;
    SharedArray<T> m_data;
    IntegralArray2e(size_t norb, size_t size): IntegralArray<T>(norb), m_data(size){}

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
    using typename IntegralArray<T>::real_t;
    using typename IntegralArray<T>::cmplx_t;
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
    using typename IntegralArray<T>::real_t;
    using typename IntegralArray<T>::cmplx_t;
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
 */

template<typename T>
struct IntegralArray2e_ih : IntegralArray2e<T> {
    using IntegralArrayBase::trig;
    using typename IntegralArray<T>::real_t;
    using typename IntegralArray<T>::cmplx_t;
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