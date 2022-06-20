//
// Created by Robert J. Anderson on 08/08/2021.
//

#ifndef M7_INTEGRALSTORAGE_H
#define M7_INTEGRALSTORAGE_H

#include <cstdlib>
#include "M7_lib/defs.h"
#include <M7_lib/parallel/SharedArray.h>

struct IntegralIndexer {
    const size_t m_norb, m_size;
    IntegralIndexer(size_t norb, size_t size): m_norb(norb), m_size(size){}
};

template<typename T, typename storage_t>
struct IntegralStorage {
    static_assert(
            std::is_same<storage_t, std::vector<T>>::value ||
            std::is_same<storage_t, SharedArray<T>>::value, "storage class must be either private or shared");
    storage_t m_data;
    const size_t m_size;
    IntegralStorage(size_t size): m_data(size), m_size(size){}

private:
    bool set_data(std::vector<T>& data, size_t iflat, T elem){
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        auto& ref = data[iflat];
        if (ref!=T(0) && !consts::nearly_equal(elem, ref, defs::helem_tol)) return false;
        ref = elem;
        return true;
    }

    bool set_data(SharedArray<T>& data, size_t iflat, T elem){
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        if (data[iflat]!=T(0) && !consts::nearly_equal(elem, data[iflat], defs::helem_tol)) return false;
        data.set(iflat, elem);
        return true;
    }

public:
    bool set_data(size_t iflat, T elem){
        return set_data(m_data, iflat, elem);
    }

    T get_data(size_t iflat) const {
        DEBUG_ASSERT_LT(iflat, m_size, "flat index OOB");
        return m_data[iflat];
    }
};

template<typename T>
using PrivateIntegralStorage = IntegralStorage<T, std::vector<T>>;
template<typename T>
using SharedIntegralStorage = IntegralStorage<T, SharedArray<T>>;

#endif //M7_INTEGRALSTORAGE_H
