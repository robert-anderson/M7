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
    MPI_Win m_win;
    T* m_data = nullptr;
public:
    SharedArray(size_t size):m_size(size) {
#ifdef HAVE_MPI
        if(mpi::on_node_i_am_root()){
            auto ierr = MPI_Win_allocate_shared(size*sizeof(T), sizeof(T), MPI_INFO_NULL, MPI_COMM_WORLD, (void*)&m_data, &m_win);
            if (ierr) throw std::runtime_error("MPI Shared memory error");
            memset(m_data, 0, size*sizeof(T));
        }
        else {
            auto ierr = MPI_Win_allocate_shared(0, sizeof(T), MPI_INFO_NULL, MPI_COMM_WORLD, (void*)&m_data, &m_win);
            if (ierr) throw std::runtime_error("MPI Shared memory error");
            int disp_unit;
            MPI_Aint alloc_size;
            MPI_Win_shared_query(m_win, 0, &alloc_size, &disp_unit, (void*)&m_data);
            ASSERT(disp_unit==sizeof(T))
            ASSERT((size_t)alloc_size==(size*sizeof(T)))
        }

#else
        m_data = new T[](size);
#endif
    }

    ~SharedArray(){
#ifdef HAVE_MPI
        MPI_Win_free(&m_win);
#else
        delete m_data;
#endif
    }

    const size_t &size(){
        return m_size;
    }

    T& operator[](const size_t &i){
        ASSERT(i<m_size)
        return *(m_data+i);
    }

    const T& operator[](const size_t &i) const{
        ASSERT(i<m_size)
        return *(m_data+i);
    }
};


#endif //M7_SHAREDARRAY_H
