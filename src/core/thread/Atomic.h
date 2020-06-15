//
// Created by Robert John Anderson on 2020-04-09.
//

#ifndef M7_ATOMIC_H
#define M7_ATOMIC_H

#include <complex>
#include "src/core/util/consts.h"

template<typename T>
struct Atomic {
    T &m_v;

    Atomic(T &v) : m_v(v) {}

    T &operator+=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
    }

    T &operator-=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
    }

    T &operator*=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
    }

    T &operator/=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
    }
};

template<>
struct Atomic<bool> {
    bool &m_v;
    Atomic(bool &v) : m_v(v) {}

    bool &operator&=(const bool &rhs) {
#pragma omp atomic update
        m_v &= rhs;
    }

    bool &operator|=(const bool &rhs) {
#pragma omp atomic update
        m_v |= rhs;
    }
};

template<typename T>
struct Atomic<std::complex<T>> {
    T &m_real;
    T &m_imag;

    Atomic(std::complex<T> &v) :
        m_real(reinterpret_cast<T(&)[2]>(v)[0]),
        m_imag(reinterpret_cast<T(&)[2]>(v)[1]){}

    std::complex<T> &operator+=(const std::complex<T> &rhs) {
        const T& real = rhs.real(); const T& imag = rhs.imag();
#pragma omp atomic update
        m_real += real;
#pragma omp atomic update
        m_imag += imag;
    }

    std::complex<T> &operator-=(const std::complex<T> &rhs) {
        const T& real = rhs.real(); const T& imag = rhs.imag();
#pragma omp atomic update
        m_real -= real;
#pragma omp atomic update
        m_imag -= imag;
    }

    std::complex<T> &operator*=(const std::complex<T> &rhs) {
        const T& real = rhs.real(); const T& imag = rhs.imag();
#pragma omp atomic update
        m_real *= real;
#pragma omp atomic update
        m_imag *= imag;
    }

    std::complex<T> &operator/=(const std::complex<T> &rhs) {
        const T& real = rhs.real(); const T& imag = rhs.imag();
#pragma omp atomic update
        m_real /= real;
#pragma omp atomic update
        m_imag /= imag;
    }

};

template<typename T>
static Atomic<T> as_atomic(T &v) {
    return Atomic<T>(v);
}

#endif //M7_ATOMIC_H
