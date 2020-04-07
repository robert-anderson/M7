//
// Created by rja on 11/03/2020.
//

#include <cstddef>
#include <array>

#ifndef M7_ARRAYINDEXER_H
#define M7_ARRAYINDEXER_H


template<size_t nind>
class ArrayIndexer {
protected:
    std::array<size_t, nind> m_shape{};
    std::array<size_t, nind> m_strides{};
public:

    explicit ArrayIndexer(const std::array<size_t, nind> &shape) {
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

    bool operator==(const ArrayIndexer &rhs) const {
        return m_shape == rhs.m_shape;
    }

    bool operator!=(const ArrayIndexer &rhs) const {
        return !(rhs == *this);
    }

private:
    template<typename T>
    void set_shape(const T &extent) {
        static_assert(std::is_integral<T>::value, "Shape requires an integral extent.");
        m_shape.fill(extent);
    }

    void set_strides() {
        m_strides.back() = 1ul;
        for (auto i = 2ul; i <= nind; i++) {
            m_strides[nind - i] = m_strides[nind - i + 1] * m_shape[nind - i + 1];
        }
    }

public:

    size_t get(const std::array<size_t, nind> &inds) const {
        size_t result = 0ul;
        for (size_t i=0ul; i<nind; ++i) result+=inds[i]*m_strides[i];
        return result;
    }
};

#endif //M7_ARRAYINDEXER_H
