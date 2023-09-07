//
// Created by Robert John Anderson on 06/09/2023.
//

#ifndef M7_LOCK_H
#define M7_LOCK_H

#include "MPIWrapper.h"

namespace lock {

    struct Vector {
        v_t<MPI_Win> m_wins;
        void resize(uint_t n) {
            for (auto& win: m_wins) MPI_Win_free(&win);
            mpi::barrier_on_node();
            m_wins.resize(n);
            uint_t buffer;
            auto ptr = reinterpret_cast<void*>(&buffer);
            for (auto& win: m_wins)
                MPI_Win_allocate(0, sizeof(int), MPI_INFO_NULL, mpi::g_node_comm, ptr, &win);
        }

        void acquire(uint_t i) {
            const auto win = m_wins[i];
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        }

        void free(uint_t i) {
            const auto win = m_wins[i];
            MPI_Win_unlock(0, win);
        }

    };
}

#endif //M7_LOCK_H
