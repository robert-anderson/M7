//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY1E_H
#define M7_INTEGRALARRAY1E_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/defs.h>
#include "IntegralArray.h"

template<typename T>
struct IntegralArray1e : IntegralArray<T> {
    using typename IntegralArray<T>::real_t;
    using typename IntegralArray<T>::cmplx_t;
    typedef std::vector<real_t> real_array_t;
    typedef std::vector<cmplx_t> cmplx_array_t;
    std::vector<T> m_data;

    IntegralArray1e(size_t norb, size_t size) : IntegralArray<T>(norb), m_data(size, 0.0) {}

protected:
    static T get_data(const std::vector<T> &data, size_t iflat) {
        DEBUG_ASSERT_LT(iflat, data.size(), "flat 1e integral index OOB");
        return data[iflat];
    }

    void set_data(std::vector<T> &data, size_t iflat, T elem) {
        DEBUG_ASSERT_LT(iflat, data.size(), "flat 1e integral index OOB");
        auto& existing = data[iflat];
        REQUIRE_TRUE(existing==0.0 || consts::nearly_equal(existing, elem),
                     "Overwriting integral data with different value: wrong permutational symmetries assumed?");
        existing = elem;
    }


private:
    virtual void ref_get(size_t i, size_t j, const real_array_t &data, real_t &elem) const = 0;

    virtual void ref_get(size_t i, size_t j, const cmplx_array_t &data, cmplx_t &elem) const = 0;

    virtual void ref_set(size_t i, size_t j, real_array_t &data, const real_t &elem) = 0;

    virtual void ref_set(size_t i, size_t j, cmplx_array_t& data, const cmplx_array_t &elem) = 0;

public:
    T get(size_t i, size_t j) const {
        defs::ham_t elem;
        ref_get(i, j, m_data, elem);
        return elem;
    }

    void set(size_t i, size_t j, const defs::ham_t &elem) {
        ref_set(i, j, m_data, elem);
    }
};

template<typename T>
struct IntegralArray1e_nosym : IntegralArray1e<T> {
    using IntegralArray1e<T>::m_norb;
    using typename IntegralArray1e<T>::real_t;
    using typename IntegralArray1e<T>::real_array_t;
    using typename IntegralArray1e<T>::cmplx_t;
    using typename IntegralArray1e<T>::cmplx_array_t;

    IntegralArray1e_nosym(size_t norb) : IntegralArray1e<T>(norb, utils::pow<2>(norb)) {}

private:
    size_t flatten(size_t i, size_t j) const {
        return i * m_norb + j;
    }
protected:

    void get(size_t i, size_t j, const real_array_t &data, real_t &elem) const override {
        elem = get_data(data, flatten(i, j));
    }

    void get(size_t i, size_t j, const cmplx_array_t &data, cmplx_t &elem) const override {
        elem = get_data(data, flatten(i, j));
    }

    void set(size_t i, size_t j, real_array_t &data, const real_t &elem) override {
        set_data(data, flatten(i, j));
    }

    void set(size_t i, size_t j, cmplx_array_t& data, const cmplx_t &elem) override {
        set_data(data, flatten(i, j));
    }
};

template<typename T>
struct IntegralArray1e_h : IntegralArray1e<T> {
    using IntegralArrayBase::trig;
    using IntegralArray1e<T>::m_norb;
    using typename IntegralArray1e<T>::real_t;
    using typename IntegralArray1e<T>::real_array_t;
    using typename IntegralArray1e<T>::cmplx_t;
    using typename IntegralArray1e<T>::cmplx_array_t;

    IntegralArray1e_h(size_t norb) : IntegralArray1e<T>(norb, IntegralArrayBase::trig(norb, 0)) {}

private:
    size_t flatten(size_t i, size_t j) const {
        return i >= j ? trig(i, j) : trig(j, i);
    }
protected:

    void get(size_t i, size_t j, const real_array_t &data, real_t &elem) const override {
        elem = get_data(data, flatten(i, j));
    }

    void get(size_t i, size_t j, const cmplx_array_t &data, cmplx_t &elem) const override {
        elem = get_data(data, flatten(i, j));
    }

    void set(size_t i, size_t j, real_array_t &data, const real_t &elem) override {
        set_data(data, flatten(i, j), elem);
    }

    void set(size_t i, size_t j, cmplx_array_t &data, const cmplx_t &elem) override {
        set_data(data, flatten(i, j), elem);
    }
};


#endif //M7_INTEGRALARRAY1E_H
