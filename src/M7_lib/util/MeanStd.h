//
// Created by rja on 23/06/22.
//

#ifndef M7_MEANSTD_H
#define M7_MEANSTD_H

#include <vector>
#include "Math.h"

template<typename T>
struct MeanStd {
    T m_mean {};
    T m_std {};

    MeanStd(typename v_t<T>::const_iterator begin, typename v_t<T>::const_iterator end) {
        // accumulate first the square mean in the m_std member
        for (auto i = begin; i != end; i++) {
            m_mean += *i;
            m_std += math::pow<2>(*i);
        }
        m_mean /= std::distance(begin, end);
        m_std /= std::distance(begin, end);
        m_std = std::sqrt(std::abs(m_std - math::pow<2>(m_mean)));
    }

    MeanStd(const v_t<T>& v): MeanStd(v.cbegin(), v.cend()){}

    static MeanStd<T> product(const MeanStd<T> &a, const MeanStd<T> &b) {
        /*
         * combine statistics in quadrature for product of random variables
         * a and b assuming no covariance
         */
        return {
            a.m_mean * b.m_mean,
            std::abs(a.m_mean * b.m_mean) * std::sqrt(math::pow<2>(a.m_std / a.m_mean) + math::pow<2>(b.m_std / b.m_mean))
        };
    }
    
    static MeanStd<T> quotient(const MeanStd<T> &a, const MeanStd<T> &b) {
        /*
         * combine statistics in quadrature for quotient of random variables
         * a and b assuming no covariance
         */
        return {
            a.m_mean / b.m_mean,
            std::abs(a.m_mean / b.m_mean) * std::sqrt(math::pow<2>(a.m_std / a.m_mean) + math::pow<2>(b.m_std / b.m_mean))
        };
    }
};


#endif //M7_MEANSTD_H
