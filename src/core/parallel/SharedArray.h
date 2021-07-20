//
// Created by rja on 13/07/2020.
//

#ifndef M7_SHAREDARRAY_H
#define M7_SHAREDARRAY_H

#include <cstddef>
#include "MPIWrapper.h"

template<typename T>
class SharedArray {
    const size_t m_size;
#ifdef ENABLE_MPI
    MPI_Win m_win;
#endif
    T *m_data = nullptr;
public:
    SharedArray(size_t size) : m_size(size) {
        /*
         * MPI_Aint window_size; double *window_data; MPI_Win node_window;
         * if (onnode_procid==0)
         * window_size = sizeof(double);
         * else window_size = 0;
         * MPI_Win_allocate_shared (window_size,sizeof(double),MPI_INFO_NULL, nodecomm, &window_data,&node_window)
         */
#ifdef ENABLE_MPI
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
        ASSERT((size_t) alloc_size == size * sizeof(T))
        MPI_Win_unlock_all(m_win);
        if (mpi::on_node_i_am_root()) std::memset(m_data, 0, size * sizeof(T));
        mpi::barrier_on_node();
#else
        m_rows = new T[m_size];
        std::memset((void*)m_rows, 0, m_size * sizeof(T));
#endif
    }

    SharedArray(SharedArray &&rhs) : m_size(rhs.m_size) {
        m_data = rhs.m_data;
#ifdef ENABLE_MPI
        m_win = rhs.m_win;
        /*
         * nullify memory window handle in rhs so that the destructor does not
         * free memory still in use here
         */
        rhs.m_win = MPI_WIN_NULL;
#else
        rhs.m_rows = nullptr;
#endif
    }

    SharedArray(const SharedArray &rhs) : SharedArray(rhs.m_size) {}

    ~SharedArray() {
        ASSERT(m_data)
#ifdef ENABLE_MPI
        if (m_win != MPI_WIN_NULL) MPI_Win_free(&m_win);
#else
        if (m_rows) delete m_rows;
#endif
    }

    const size_t &size() {
        return m_size;
    }

    void set(const size_t &i, const T &v) {
        // element-modifying access should only take place on the root rank
        if (mpi::on_node_i_am_root()) {
            ASSERT(i < m_size)
            m_data[i] = v;
        }
    }

    const T &get(const size_t &i) const {
        ASSERT(i < m_size)
        return m_data[i];
    }

    const T &operator[](const size_t &i) const {
        return get(i);
    }
};


#endif //M7_SHAREDARRAY_H
