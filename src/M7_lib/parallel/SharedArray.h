//
// Created by Robert J. Anderson on 13/07/2020.
//

#ifndef M7_SHAREDARRAY_H
#define M7_SHAREDARRAY_H

#include <cstddef>
#include "MPIWrapper.h"
#include "MPIAssert.h"

template<typename T>
class SharedArray {
    const uint_t m_size;
    MPI_Win m_win;
    T *m_data = nullptr;
public:
    SharedArray(uint_t size) : m_size(size) {
        auto baseptr = reinterpret_cast<void *>(&m_data);
        if (mpi::on_node_i_am_root()) {
            auto ierr = MPI_Win_allocate_shared(size * sizeof(T), sizeof(T), MPI_INFO_NULL, g_node_comm,
                                                baseptr, &m_win);
            if (ierr) throw std::runtime_error("MPI Shared memory error");
        } else {
            auto ierr = MPI_Win_allocate_shared(0, sizeof(T), MPI_INFO_NULL, g_node_comm, baseptr, &m_win);
            if (ierr) throw std::runtime_error("MPI Shared memory error");
        }
        MPI_Win_lock_all(0, m_win);
        MPI_Win_sync(m_win);
        mpi::barrier_on_node();
        int disp_unit;
        MPI_Aint alloc_size;
        /*
         * MPI_Aint window_size0; int window_unit; double *win0_addr;
         * MPI_Win_shared_query(node_window, 0, &window_size0, &window_unit, &win0_addr);
         */
        auto ierr = MPI_Win_shared_query(m_win, 0, &alloc_size, &disp_unit, baseptr);
        if (ierr != MPI_SUCCESS) throw std::runtime_error("MPI Memory Window query failed");
        ASSERT(disp_unit == sizeof(T))
        ASSERT((uint_t) alloc_size == size * sizeof(T))
        MPI_Win_unlock_all(m_win);
        if (mpi::on_node_i_am_root()) std::memset(reinterpret_cast<char*>(m_data), 0, size * sizeof(T));
        mpi::barrier_on_node();
    }

    SharedArray(SharedArray &&rhs) : m_size(rhs.m_size) {
        m_data = rhs.m_data;
        m_win = rhs.m_win;
        /*
         * nullify memory window handle in rhs so that the destructor does not
         * free memory still in use here
         */
        rhs.m_win = MPI_WIN_NULL;
    }

    SharedArray(const SharedArray &rhs) : SharedArray(rhs.m_size) {}

    ~SharedArray() {
        ASSERT(m_data)
        if (m_win != MPI_WIN_NULL) MPI_Win_free(&m_win);
    }

    const uint_t &size() {
        return m_size;
    }

    void set(const uint_t &i, const T &v) {
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) {
            DEBUG_ASSERT_LT(i, m_size, "SharedArray element OOB");
            m_data[i] = v;
        }
    }

    const T &operator[](const uint_t &i) const {
        DEBUG_ASSERT_LT(i, m_size, "SharedArray element OOB");
        return m_data[i];
    }
};


#endif //M7_SHAREDARRAY_H
