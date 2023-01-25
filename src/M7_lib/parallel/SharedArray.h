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
    uint_t m_nelement = 0;
    const uint_t m_element_size;
    uint_t m_nbyte = 0;
    void *m_data = nullptr;
private:

    SharedArrayBase(uint_t element_size);

    static void alloc(uint_t nelement, uint_t element_size, MPI_Win* win, void** data);

    static void free(MPI_Win* win, void** data);

    void alloc(uint_t nelement, uint_t element_size);

    void free();



protected:
    MPI_Win m_win;
public:
    SharedArrayBase(uint_t nelement, uint_t element_size);

    SharedArrayBase& operator=(const SharedArrayBase& other);

    SharedArrayBase& operator=(SharedArrayBase&& other);

    SharedArrayBase(const SharedArrayBase &other);

    SharedArrayBase(SharedArrayBase &&other);

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
