//
// Created by rja on 19/05/2020.
//

#ifndef M7_REDUCTION_H
#define M7_REDUCTION_H

//enum

#include <limits>

#if 0
namespace reduction {
    template<typename T>
    void max(const T& priv, T& shared){
        shared = std::numeric_limits<T>::lowest();
#pragma omp parallel critical (reduce_max)
        {
            if (priv>shared) shared=priv;
        }
    }
    template<typename T>
    void min(const T& priv, T& shared){
        shared = std::numeric_limits<T>::max();
#pragma omp parallel critical (reduce_min)
        {
            if (priv<shared) shared=priv;
        }
    }
}
#endif


#endif //M7_REDUCTION_H
