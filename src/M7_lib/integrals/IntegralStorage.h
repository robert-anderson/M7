//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALSTORAGE_H
#define M7_INTEGRALSTORAGE_H

#include <cstdlib>
#include "M7_lib/defs.h"
#include <M7_lib/parallel/SharedArray.h>

struct IntegralIndexer {
    const uint_t m_norb, m_size;
    IntegralIndexer(uint_t norb, uint_t size): m_norb(norb), m_size(size){}
};

template<typename T, typename storage_t>
struct IntegralStorage {
    static_assert(
            std::is_same<storage_t, v_t<T>>::value ||
            std::is_same<storage_t, SharedArray<T>>::value, "storage class must be either private or shared");
    storage_t m_data;
    const uint_t m_size;
    static constexpr arith::comp_t<T> c_coeff_atol = 1e-9;

    IntegralStorage(uint_t size): m_data(size), m_size(size){}

private:
    /**
     * @param elem
     *  coefficient element
     * @return
     *  true if the magnitude of elem exceeds minimum
     */
    static bool is_significant(T elem) {
        return std::abs(elem) >= c_coeff_atol;
    }

    bool set_data_(v_t<T>& data, uint_t iflat, T elem){
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        if (!is_significant(elem)) return true;
        auto& ref = data[iflat];
        if (ref!=T(0) && is_significant(elem-ref)) return false;
        ref = elem;
        return true;
    }

    bool set_data_(SharedArray<T>& data, uint_t iflat, T elem){
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        if (!is_significant(elem)) return true;
        if (data[iflat]!=T(0) && is_significant(elem-data[iflat])) return false;
        data.set_(iflat, elem);
        return true;
    }

public:
    bool set_data_(uint_t iflat, T elem){
        DEBUG_ASSERT_TRUE(mpi::on_node_i_am_root(), "integral elements should only be set from node root ranks");
        return set_data_(m_data, iflat, elem);
    }

    T get_data(uint_t iflat) const {
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        return m_data[iflat];
    }
};

/*
 * private storage is unused, all integral data structures are shared over the node (machine or NUMA region)
 */
//template<typename T>
//using PrivateIntegralStorage = IntegralStorage<T, v_t<T>>;
template<typename T>
using SharedIntegralStorage = IntegralStorage<T, SharedArray<T>>;

#endif //M7_INTEGRALSTORAGE_H
