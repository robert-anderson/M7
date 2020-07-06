//
// Created by Robert John Anderson on 2020-04-09.
//

#ifndef M7_ATOMIC_H
#define M7_ATOMIC_H

#include <complex>
#include "src/core/util/consts.h"
#include "omp.h"

template<typename T>
struct Atomic {
    T &m_v;

    Atomic(T &v) : m_v(v) {}

    T &operator+=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
        return m_v;
    }

    T &operator-=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
        return m_v;
    }

    T &operator*=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
        return m_v;
    }

    T &operator/=(const T &rhs) {
#pragma omp atomic update
        m_v += rhs;
        return m_v;
    }
};

template<>
struct Atomic<bool> {
    bool &m_v;
    Atomic(bool &v) : m_v(v) {}

    bool &operator&=(const bool &rhs) {
#pragma omp atomic update
        m_v &= rhs;
        return m_v;
    }

    bool &operator|=(const bool &rhs) {
#pragma omp atomic update
        m_v |= rhs;
        return m_v;
    }
};

template<typename T>
struct Atomic<std::complex<T>> {
    std::complex<T> &m_v;

    Atomic(std::complex<T> &v) : m_v(v){}

    std::complex<T> &operator+=(const std::complex<T> &rhs) {
        T& real = reinterpret_cast<T(&)[2]>(m_v)[0];
        T& imag = reinterpret_cast<T(&)[2]>(m_v)[1];
        const T& delta_real = rhs.real(); const T& delta_imag = rhs.imag();
#pragma omp atomic update
        real += delta_real;
#pragma omp atomic update
        imag += delta_imag;
        return m_v;
    }

    std::complex<T> &operator-=(const std::complex<T> &rhs) {
        T& real = reinterpret_cast<T(&)[2]>(m_v)[0];
        T& imag = reinterpret_cast<T(&)[2]>(m_v)[1];
        const T& delta_real = rhs.real(); const T& delta_imag = rhs.imag();
#pragma omp atomic update
        real -= delta_real;
#pragma omp atomic update
        imag -= delta_imag;
        return m_v;
    }

    std::complex<T> &operator*=(const std::complex<T> &rhs) {
        T& real = reinterpret_cast<T(&)[2]>(m_v)[0];
        T& imag = reinterpret_cast<T(&)[2]>(m_v)[1];
        const T& delta_real = rhs.real(); const T& delta_imag = rhs.imag();
#pragma omp atomic update
        real *= delta_real;
#pragma omp atomic update
        imag *= delta_imag;
        return m_v;
    }

    std::complex<T> &operator/=(const std::complex<T> &rhs) {
        T& real = reinterpret_cast<T(&)[2]>(m_v)[0];
        T& imag = reinterpret_cast<T(&)[2]>(m_v)[1];
        const T& delta_real = rhs.real(); const T& delta_imag = rhs.imag();
#pragma omp atomic update
        real /= delta_real;
#pragma omp atomic update
        imag /= delta_imag;
        return m_v;
    }

};

template<typename T>
static Atomic<T> as_atomic(T &v) {
    return Atomic<T>(v);
}

#endif //M7_ATOMIC_H
