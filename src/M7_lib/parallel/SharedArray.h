//
// Created by Robert J. Anderson on 13/07/2020.
//

#ifndef M7_SHAREDARRAY_H
#define M7_SHAREDARRAY_H

#include <cstddef>
#include "MPIWrapper.h"
#include "MPIAssert.h"

class SharedArrayBase {
public:
    const uint_t m_nelement;
    const uint_t m_element_size;
    const uint_t m_nbyte;
    void *m_data = nullptr;
protected:
    MPI_Win m_win;
public:
    SharedArrayBase(uint_t nelement, uint_t element_size);

    SharedArrayBase(SharedArrayBase &&rhs);

    SharedArrayBase(const SharedArrayBase &rhs) : SharedArrayBase(rhs.m_nelement, rhs.m_element_size) {}

    ~SharedArrayBase();
};

template<typename T>
class SharedArray : public SharedArrayBase {
public:
    SharedArray(uint_t size) : SharedArrayBase(size, sizeof(T)) {}

    uint_t size() const {
        return m_nelement;
    }

    void set(uint_t i, const T &v) {
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) {
            DEBUG_ASSERT_LT(i, size(), "SharedArray element OOB");
            reinterpret_cast<T*>(m_data)[i] = v;
        }
    }

    const T &operator[](uint_t i) const {
        DEBUG_ASSERT_LT(i, size(), "SharedArray element OOB");
        return reinterpret_cast<const T*>(m_data)[i];
    }
};


#endif //M7_SHAREDARRAY_H
