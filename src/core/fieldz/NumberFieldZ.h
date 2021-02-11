//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELDZ_H
#define M7_NUMBERFIELDZ_H

#include "FieldBaseZ.h"
#include "src/core/nd/NdFormat.h"

template<typename T, size_t nind>
struct NumberFieldZ : FieldBaseZ {
    NdFormat<nind> m_format;

    NumberFieldZ(std::array<size_t, nind> shape) :
            FieldBaseZ(sizeof(T) * NdFormat<nind>(shape).nelement(), typeid(NumberFieldZ<T, nind>)),
            m_format(shape) {}

    NumberFieldZ &operator=(const T &v) {
        std::fill((T*)raw_view(), ((T*)raw_view())+m_element_size, v);
        return *this;
    }

    NumberFieldZ &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == m_format.nelement());
        std::memcpy(raw_view(), v.data(), m_element_size);
        return *this;
    }

    const T &operator[](const size_t &i) const {
        ASSERT(i < m_format.nelement());
        return ((T *) raw_view())[i];
    }

    T &operator[](const size_t &i) {
        ASSERT(i < m_format.nelement());
        return ((T *) raw_view())[i];
    }

    template<typename ...Args>
    T &operator()(Args... inds) {
        return ((T *) raw_view())[m_format.flatten(inds...)];
    }

    template<typename ...Args>
    const T &operator()(Args... inds) const {
        return ((T *) raw_view())[m_format.flatten(inds...)];
    }

    std::string to_string() const override {
        ASSERT(m_view_offset!=~0ul)
        std::string tmp = "[";
        const auto nelement = m_format.nelement();
        for (size_t i = 0ul; i < nelement; ++i)
            tmp += std::to_string(this->operator[](i)) + " ";
        return tmp + "]";
    }
};


struct BosonOnvFieldZ : NumberFieldZ<uint8_t, 1> {
    size_t m_nmode;

    BosonOnvFieldZ(size_t nmode) :
            NumberFieldZ<uint8_t, 1>({nmode}), m_nmode(nmode) {}
};



#endif //M7_NUMBERFIELDZ_H
