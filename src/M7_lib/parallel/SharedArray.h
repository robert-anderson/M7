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
    /**
     * The rank index in the global communicator that has the exclusive right to write on m_data
     */
    const uint_t m_irank_owner;
private:

    static void alloc(uint_t nelement, uint_t element_size, MPI_Win* win, void** data);

    static void free(MPI_Win* win, void** data);

    void alloc(uint_t nelement);

    void free();

protected:
    MPI_Win m_win;

    SharedArrayBase(uint_t element_size, uint_t irank_owner);

    // owner defaults to the root rank of the shared memory region
    explicit SharedArrayBase(uint_t element_size);
public:

    SharedArrayBase(): SharedArrayBase(1ul){}

    SharedArrayBase(uint_t nelement, uint_t element_size, uint_t irank_owner);

    SharedArrayBase& operator=(const SharedArrayBase& other);

    SharedArrayBase& operator=(SharedArrayBase&& other);

    SharedArrayBase(const SharedArrayBase &other);

    SharedArrayBase(SharedArrayBase &&other);

    ~SharedArrayBase();

protected:
    void set_(uint_t i, uint_t n, const void* src) {
        DEBUG_ASSERT_TRUE(mpi::i_am_root(mpi::SharedMemory), "element-modifying access should only take place on the root rank");
        DEBUG_ASSERT_LT(i, m_nelement, "begin OOB");
        DEBUG_ASSERT_LE(i+n, m_nelement, "end OOB");
        std::memcpy(m_data+(i*m_element_size), src, n*m_element_size);

    }

    void set_(uint_t i, const void* src) {
        set_(i, 1, src);
    }

    void set_(const void* src) {
        set_(0, m_nelement, src);
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
    SharedArray(uint_t size, uint_t irank_owner) : SharedArrayBase(size, sizeof(T), irank_owner) {}
    SharedArray(uint_t size) : SharedArray(size, mpi::irank_world_shmem_root()){}

    uint_t size() const {
        return m_nelement;
    }

    void set_(uint_t i, const T &v) {
        SharedArrayBase::set_(i, &v);
    }

    void set_(uint_t i, const v_t<T> &v) {
        SharedArrayBase::set_(i, v.size(), v.data());
    }

    void set_(const v_t<T> &v) {
        SharedArrayBase::set_(v.data());
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

    const T* cbegin() const {
        return reinterpret_cast<const T*>(m_data);
    }
};


template<typename T>
class SharedScalar : protected SharedArray<T> {

public:
    SharedScalar(uint_t irank_owner) : SharedArray<T>(1, irank_owner){}
    SharedScalar() : SharedScalar(mpi::irank_world_shmem_root()){}

    explicit SharedScalar(uint_t irank_owner, const T& v) : SharedScalar(irank_owner) {
        if (mpi::i_am_root(mpi::SharedMemory)) set_(v);
        mpi::barrier(mpi::SharedMemory);
    }

    void set_(const T &v) {
        SharedArrayBase::set_(0, &v);
    }

    operator const T& () const {
        return reinterpret_cast<const T&>(*SharedArrayBase::m_data);
    }
};

#endif //M7_SHAREDARRAY_H
