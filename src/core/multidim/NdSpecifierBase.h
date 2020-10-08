/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 *
 * A MultiDimensional
 *
 */



#ifndef M7_NDSPECIFIERBASE_H
#define M7_NDSPECIFIERBASE_H

#include "cstddef"
#include "array"

template<typename selector_t, size_t nind>
class NdSpecifierBase {
    selector_t m_selector;
protected:
    std::array<size_t, nind> m_shape{};
    std::array<size_t, nind> m_strides{};
public:

    template<typename ...Args>
    NdSpecifierBase(selector_t selector, const size_t &first, Args ...shape):
            m_selector(selector) {
        static_assert(sizeof...(shape) + 1 == nind, "Invalid number of shape arguments.");
        set_shape(first, shape...);
        set_strides();
    }

    NdSpecifierBase(const std::array<size_t, nind> &shape) {
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

public:
    typename selector_t::accessor_t select(const size_t &flat) {
        return m_selector(flat);
    }

    typename selector_t::const_accessor_t select(const size_t &flat) const {
        return m_selector(flat);
    }

    std::string to_string() const {
        std::string out;
        for (size_t i = 0ul; i < nelement(); ++i) {
            out += " " + select(i).to_string();
        }
        return out;
    }
};


#endif //M7_NDSPECIFIERBASE_H
