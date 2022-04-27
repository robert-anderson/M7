//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY2E_H
#define M7_INTEGRALARRAY2E_H

#include <M7_lib/defs.h>
#include <M7_lib/parallel/SharedArray.h>

#include "IntegralArray.h"

/**
 * store two-electron integrals with the chemists' notation (ij|kl) where the integrand is i*(x1)j(x1) 1/r12 k*(x2)l(x2)
 * such quantities admit varying levels of permutational symmetry whereby an integral indexed by (ij|kl) is equal (up to
 * a complex conjugation) to any of the following
 *
 *  (kl|ji)     (fermion indistinguishability "I")
 *  (ji|lk)     (hermiticity "H")
 *  (ij|lk)     (real orbitals "R")
 *
 * The assumption of:
 *  - none of the above is called "1-fold" symmetry,
 *  - of "I" only is called "2-fold" symmetry,
 *  - of "IH" is called "4 fold" symmetry,
 *  - and of "IHR" is called "8 fold symmetry"
 *
 *  1 fold makes no saving in terms of memory overhead, and thus places no constraints on the indices stored
 *  2 fold computes combined ij and kl indices, and only stores the elements with ij>=kl
 *  8 fold only stores the elements with i>=j and k>=l and ij>=kl
 *  4 fold is a slightly more complicated case, which is important in situations where the integrals correspond to
 *      complex-valued orbitals. the full explanation is given in the subclass docs
 */
struct IntegralArray2e : IntegralArray {
    SharedArray<defs::ham_t> m_data;
    IntegralArray2e(size_t norb, size_t size):
        IntegralArray(norb), m_data(size){}

private:
    virtual void get(const size_t& i, const size_t& j, const size_t& k, const size_t& l,
                     const SharedArray<defs::ham_comp_t> &data, defs::ham_comp_t& elem) const = 0;
    virtual void get(const size_t& i, const size_t& j, const size_t& k, const size_t& l,
                     const SharedArray<std::complex<defs::ham_comp_t>> &data, std::complex<defs::ham_comp_t>& elem) const = 0;
    virtual void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l,
                     SharedArray<defs::ham_comp_t> &data, const defs::ham_comp_t& elem) const = 0;
    virtual void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l,
                     SharedArray<std::complex<defs::ham_comp_t>> &data, const std::complex<defs::ham_comp_t>& elem) const = 0;

public:
    defs::ham_t get(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const {
        defs::ham_t elem;
        get(i, j, k, l, m_data, elem);
        return elem;
    }

    void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l, defs::ham_t elem) {
        set(i, j, k, l, m_data, elem);
    }
};

struct IntegralArray2e_1fold : IntegralArray2e {
    const size_t m_norb2, m_norb_3;
    IntegralArray2e_1fold(size_t norb) : IntegralArray2e(norb, utils::pow<4>(norb)),
                                         m_norb2(norb*norb), m_norb_3(m_norb2*norb){}

private:
    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        elem = data[i*m_norb_3+j*m_norb2+k*m_norb+l];
    }

    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<std::complex<defs::ham_comp_t>> &data, std::complex<defs::ham_comp_t> &elem) const override {
        elem = data[i*m_norb_3+j*m_norb2+k*m_norb+l];
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, SharedArray<defs::ham_comp_t> &data,
             const defs::ham_comp_t &elem) const override {
        data.set(i*m_norb_3+j*m_norb2+k*m_norb+l, elem);
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             SharedArray<std::complex<defs::ham_comp_t>> &data,
             const std::complex<defs::ham_comp_t> &elem) const override {
        data.set(i*m_norb_3+j*m_norb2+k*m_norb+l, elem);
    }
};

struct IntegralArray2e_2fold : IntegralArray2e {
    IntegralArray2e_2fold(size_t norb) : IntegralArray2e(norb, trig(norb*norb, 0)){}

private:
    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        auto ij = i*m_norb+j;
        auto kl = k*m_norb+l;
        elem = ij>=kl ? data[trig(ij, kl)] : data[trig(kl, ij)];
    }

    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<std::complex<defs::ham_comp_t>> &data, std::complex<defs::ham_comp_t> &elem) const override {
        auto ij = i*m_norb+j;
        auto kl = k*m_norb+l;
        elem = ij>=kl ? data[trig(ij, kl)] : data[trig(kl, ij)];
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, SharedArray<defs::ham_comp_t> &data,
             const defs::ham_comp_t &elem) const override {
        auto ij = i*m_norb+j;
        auto kl = k*m_norb+l;
        data.set(ij>=kl ? trig(ij, kl) : trig(kl, ij), elem);
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             SharedArray<std::complex<defs::ham_comp_t>> &data,
             const std::complex<defs::ham_comp_t> &elem) const override {
        auto ij = i*m_norb+j;
        auto kl = k*m_norb+l;
        data.set(ij>=kl ? trig(ij, kl) : trig(kl, ij), elem);
    }
};


struct IntegralArray2e_8fold : IntegralArray2e {
    IntegralArray2e_8fold(size_t norb): IntegralArray2e(norb, trig(trig(norb, 0), trig(norb, 0))){}

private:
    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)               (lk|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ji|kl)               (kl|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk)               (lk|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
        }
    }

    /**
     * it's impossible for real orbitals to yield integrals with non-zero imaginary part, so in this case it must be
     * assumed that we are (wastefully) using a complex container for real-valued integrals
     */
    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<std::complex<defs::ham_comp_t>> &data, std::complex<defs::ham_comp_t> &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)               (lk|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ji|kl)               (kl|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk)               (lk|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
        }
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, SharedArray<defs::ham_comp_t> &data,
             const defs::ham_comp_t &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)               (lk|ij)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ji|kl)               (kl|ji)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk)               (lk|ji)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
        }
    }

    void set(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             SharedArray<std::complex<defs::ham_comp_t>> &data,
             const std::complex<defs::ham_comp_t> &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)               (lk|ij)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ji|kl)               (kl|ji)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk)               (lk|ji)
                data.set(ij>=kl ? trig(ij, kl): trig(kl, ij), elem);
            }
        }
    }
};

/**
 * 4-fold permutational symmetry is treated as though we have two tandem 8-fold arrays,
 *                                                    "I"     "H"          "IH"
 * one for the indices equivalent to (ij|kl), i.e. (kl|ij), (ji|lk), and (lk|ji)
 * and another for those equivalent to (ij|lk) i.e. (lk|ij), (ji|kl), and (kl|ji)
 */
struct IntegralArray2e_4fold : IntegralArray2e {
    const size_t m_n8fold;
    IntegralArray2e_4fold(size_t norb):
        IntegralArray2e(norb, 2*trig(trig(norb, 0), trig(norb, 0))),
        m_n8fold(m_data.size()/2){}

private:
    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)                        (lk|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)+m_n8fold]: data[trig(kl, ij)+m_n8fold];
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ji|kl)                        (kl|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)+m_n8fold]: data[trig(kl, ij)+m_n8fold];
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk)               (lk|ji)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
        }
    }

    void get(const size_t &i, const size_t &j, const size_t &k, const size_t &l,
             const SharedArray<std::complex<defs::ham_comp_t>> &data, std::complex<defs::ham_comp_t> &elem) const override {
        if (i>=j){
            auto ij = trig(i, j);
            if (k>=l){
                auto kl = trig(k, l);
                //                      (ij|kl)               (kl|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)]: data[trig(kl, ij)];
            }
            else {
                auto kl = trig(l, k);
                //                      (ij|lk)                        (lk|ij)
                elem = (ij>=kl) ? data[trig(ij, kl)+m_n8fold]: data[trig(kl, ij)+m_n8fold];
            }
        }
        else {
            auto ij = trig(j, i);
            if (k>=l){
                auto kl = trig(k, l);
                //                     (ji|kl) = (ij|lk)*                     (kl|ji) = (lk|ij)* = (ij|lk)*
                elem = std::conj((ij>=kl) ? data[trig(ij, kl)+m_n8fold]: data[trig(kl, ij)+m_n8fold]);
            }
            else {
                auto kl = trig(l, k);
                //                      (ji|lk) = (ij|kl)*               (lk|ji) = (kl|ij)* = (ij|kl)*
                elem =  std::conj((ij>=kl) ? data[trig(ij, kl)]: std::conj(data[trig(kl, ij)]));
            }
        }
    }
};





#endif //M7_INTEGRALARRAY2E_H
