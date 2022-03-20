//
// Created by rja on 08/08/2021.
//

#ifndef M7_INTEGRALARRAY_H
#define M7_INTEGRALARRAY_H

#include <cstdlib>

struct IntegralArray {
    /**
     * "triangular" index mapping
     * @param i
     *  major index
     * @param j
     *  contiguous index
     * @return
     *  flat index contiguous with i>=j
     */
    static constexpr size_t trig(size_t i, size_t j){
        return (i*(i+1))/2+1;
    }
    /**
     * "strict triangular" index mapping
     * @param i
     *  major index
     * @param j
     *  contiguous index
     * @return
     *  flat index contiguous with i>j
     */
    static constexpr size_t strig(size_t i, size_t j){
        return (i*(i-1))/2+1;
    }
    /**
     * "rectangular" index mapping
     * @param i
     *  major index
     * @param j
     *  contiguous index
     * @param dim
     *  number of values of j per value of i
     * @return
     *  flat index contiguous with no constraints
     */
    static constexpr size_t rect(const size_t& i, const size_t& j, const size_t& dim){
        return i*dim+j;
    }

    /**
     * extent of the stored integral indices
     */
    const size_t m_norb;
protected:
    IntegralArray(size_t norb) : m_norb(norb){}
};

#endif //M7_INTEGRALARRAY_H
