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
    buf_t *m_data = nullptr;
private:

    static void alloc(uint_t nelement, uint_t element_size, MPI_Win* win, void** data);

    static void free(MPI_Win* win, void** data);

    void alloc(uint_t nelement);

    void free();

protected:
    MPI_Win m_win;
public:

    SharedArrayBase(uint_t element_size);

    SharedArrayBase(): SharedArrayBase(1ul){}

    SharedArrayBase(uint_t nelement, uint_t element_size);

    SharedArrayBase& operator=(const SharedArrayBase& other);

    SharedArrayBase& operator=(SharedArrayBase&& other);

    SharedArrayBase(const SharedArrayBase &other);

    SharedArrayBase(SharedArrayBase &&other);

    ~SharedArrayBase();

protected:
    void set(uint_t i, uint_t n, const void* src) {
        DEBUG_ASSERT_LT(i, m_nelement, "begin OOB");
        DEBUG_ASSERT_LE(i+n, m_nelement, "end OOB");
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) std::memcpy(m_data+(i*m_element_size), src, n*m_element_size);

    }

    void set(uint_t i, const void* src) {
        DEBUG_ASSERT_LT(i, m_nelement, "begin OOB");
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) std::memcpy(m_data+(i*m_element_size), src, m_element_size);
    }

    void set(const void* src) {
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) std::memcpy(m_data, src, m_nbyte);
    }

    void get(uint_t i, uint_t n, void* dst) {
        DEBUG_ASSERT_LT(i, m_nelement, "begin OOB");
        DEBUG_ASSERT_LE(i+n, m_nelement, "end OOB");
        std::memcpy(dst, m_data+(i*m_element_size), n*m_element_size);

    }

    void get(uint_t i, void* dst) {
        DEBUG_ASSERT_LT(i, m_nelement, "begin OOB");
        std::memcpy(dst, m_data+(i*m_element_size), m_element_size);
    }

    void get(void* dst) {
        std::memcpy(dst, m_data, m_nbyte);
    }
};

template<typename T>
class SharedArray : public SharedArrayBase {
public:
    SharedArray(uint_t size) : SharedArrayBase(size, sizeof(T)) {}

    uint_t size() const {
        return m_nelement;
    }

    void set(uint_t i, const T &v) {
        SharedArrayBase::set(i, &v);
    }

    void set(uint_t i, const v_t<T> &v) {
        SharedArrayBase::set(i, v.size(), v.data());
    }

    void set(const v_t<T> &v) {
        DEBUG_ASSERT_TRUE(!mpi::on_node_i_am_root() || (v.size()==m_nelement), "incompatible source");
        SharedArrayBase::set(v.data());
    }

    void get(uint_t i, v_t<T> &v) {
        SharedArrayBase::get(i, v.size(), v.data());
    }

    void get(v_t<T> &v) {
        if (v.size() < m_nelement) v.resize(m_nelement);
        SharedArrayBase::get(v.data());
    }

    const T &operator[](uint_t i) const {
        DEBUG_ASSERT_LT(i, size(), "SharedArray element OOB");
        return reinterpret_cast<const T*>(m_data)[i];
    }
};


template<typename T>
class SharedScalar : protected SharedArray<T> {
    /**
     * written to if this MPI rank is not the root
     */
    T m_trash;
public:
    SharedScalar() : SharedArray<T>(1ul){}

    SharedScalar(const T& v) : SharedScalar() {
        *this = v;
    }

    SharedScalar& operator=(const T& v){
        SharedArray<T>::set(0, v);
        return *this;
    }

    operator const T& () const {
        return reinterpret_cast<const T&>(*SharedArrayBase::m_data);
    }
};

#endif //M7_SHAREDARRAY_H
