//
// Created by rja on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY1E_H
#define M7_INTEGRALARRAY1E_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/defs.h>
#include "IntegralArray.h"

struct IntegralArray1e : IntegralArray {
    std::vector<defs::ham_t> m_data;

    IntegralArray1e(size_t norb, size_t size) :
            IntegralArray(norb), m_data(size, 0.0) {}

protected:
    template<typename T>
    static const T &get_data(const std::vector<T> &data, const size_t &iflat) {
        DEBUG_ASSERT_LT(iflat, data.size(), "flat 1e integral index OOB");
        return data[iflat];
    }
    template<typename T>
    void set_data(std::vector<T> &data, const size_t &iflat, const T& elem) {
        DEBUG_ASSERT_LT(iflat, data.size(), "flat 1e integral index OOB");
        auto& existing = data[iflat];
        REQUIRE_TRUE(existing==0.0 || consts::nearly_equal(existing, elem),
                     "Overwriting integral data with different value: wrong permutational symmetries assumed?");
        existing = elem;
    }


private:
    virtual void get(const size_t &i, const size_t &j,
                     const std::vector<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const = 0;

    virtual void get(const size_t &i, const size_t &j,
                     const std::vector<std::complex<defs::ham_comp_t>> &data,
                     std::complex<defs::ham_comp_t> &elem) const = 0;

    virtual void set(const size_t &i, const size_t &j,
                     std::vector<defs::ham_comp_t> &data, const defs::ham_comp_t &elem) = 0;

    virtual void set(const size_t &i, const size_t &j,
                     std::vector<std::complex<defs::ham_comp_t>> &data,
                     const std::complex<defs::ham_comp_t> &elem) = 0;

public:
    defs::ham_t get(const size_t &i, const size_t &j) const {
        defs::ham_t elem;
        get(i, j, m_data, elem);
        return elem;
    }

    void set(const size_t &i, const size_t &j, const defs::ham_t &elem) {
        set(i, j, m_data, elem);
    }
};

struct IntegralArray1e_1fold : IntegralArray1e {
    IntegralArray1e_1fold(size_t norb) : IntegralArray1e(norb, utils::pow<2>(norb)) {}

    void get(const size_t &i, const size_t &j,
             const std::vector<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        elem = get_data(data, i * m_norb + j);
    }

    void get(const size_t &i, const size_t &j,
             const std::vector<std::complex<defs::ham_comp_t>> &data,
             std::complex<defs::ham_comp_t> &elem) const override {
        elem = get_data(data, i * m_norb + j);
    }

    void set(const size_t &i, const size_t &j, std::vector<defs::ham_comp_t> &data,
             const defs::ham_comp_t &elem) override {
        set_data(data, i * m_norb + j, elem);
    }

    void set(const size_t &i, const size_t &j, std::vector<std::complex<defs::ham_comp_t>> &data,
             const std::complex<defs::ham_comp_t> &elem) override {
        set_data(data, i * m_norb + j, elem);
    }
};

struct IntegralArray1e_2fold : IntegralArray1e {
    IntegralArray1e_2fold(size_t norb) : IntegralArray1e(norb, IntegralArray::trig(norb, 0)) {}

    void get(const size_t &i, const size_t &j,
             const std::vector<defs::ham_comp_t> &data, defs::ham_comp_t &elem) const override {
        elem = get_data(data, i >= j ? trig(i, j) : trig(j, i));
    }

    void get(const size_t &i, const size_t &j,
             const std::vector<std::complex<defs::ham_comp_t>> &data,
             std::complex<defs::ham_comp_t> &elem) const override {
        elem = (i >= j) ? get_data(data, trig(i, j)) : std::conj(get_data(data, trig(j, i)));
    }

    void set(const size_t &i, const size_t &j, std::vector<defs::ham_comp_t> &data,
             const defs::ham_comp_t &elem) override {
        i>=j ? set_data(data, trig(i, j), elem): set_data(data, trig(j, i), elem);
    }

    void set(const size_t &i, const size_t &j, std::vector<std::complex<defs::ham_comp_t>> &data,
             const std::complex<defs::ham_comp_t> &elem) override {
        i>=j ? set_data(data, trig(i, j), elem): set_data(data, trig(j, i), std::conj(elem));
    }
};


#endif //M7_INTEGRALARRAY1E_H
