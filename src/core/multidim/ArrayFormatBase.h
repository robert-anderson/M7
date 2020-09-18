//
// Created by RJA on 18/09/2020.
//

#ifndef M7_ARRAYFORMATBASE_H
#define M7_ARRAYFORMATBASE_H

#include "cstddef"
#include "array"

template<size_t nind>
class ArrayFormatBase {
protected:
    std::array<size_t, nind> m_shape{};
    std::array<size_t, nind> m_strides{};
public:

    template<typename ...Args>
    ArrayFormatBase(const size_t &first, Args ...shape) {
        static_assert(sizeof...(shape) + 1 == nind, "Invalid number of shape arguments.");
        set_shape(first, shape...);
        set_strides();
    }

    ArrayFormatBase(const std::array<size_t, nind> &shape) {
        m_shape = shape;
        set_strides();
    }

    const std::array<size_t, nind> &shape() const {
        return m_shape;
    }

    const std::array<size_t, nind> &strides() const {
        return m_strides;
    }

    size_t nelement() const {
        return m_strides.front() * m_shape.front();
    }

private:
    template<typename T>
    void set_shape(const T &extent) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape.fill(extent);
    }

    template<typename T, typename ...Args>
    void set_shape(const T &first, Args ...args) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind - sizeof...(args) - 1] = first;
        set_shape(args...);
    }

    template<typename T>
    void set_shape(const T &first, T &second) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape[nind - 2] = first;
        m_shape[nind - 1] = second;
    }

    void set_strides() {
        m_strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            m_strides[nind - i] = m_strides[nind - i + 1] * m_shape[nind - i + 1];
        }
    }
};



#endif //M7_ARRAYFORMATBASE_H
