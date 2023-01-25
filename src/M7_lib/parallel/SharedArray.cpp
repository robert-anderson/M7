//
// Created by rja on 24/01/23.
//

#include "SharedArray.h"

SharedArrayBase::SharedArrayBase(uint_t nelement, uint_t element_size) :
        m_nelement(nelement), m_element_size(element_size), m_nbyte(m_nelement * m_element_size) {
    const auto local_nbyte = mpi::on_node_i_am_root() ? m_nbyte : 0ul;
    auto ierr = MPI_Win_allocate_shared(local_nbyte, m_element_size, MPI_INFO_NULL, g_node_comm, &m_data, &m_win);
    REQUIRE_EQ(ierr, MPI_SUCCESS, "MPI Shared memory error");
    MPI_Win_lock_all(0, m_win);
    MPI_Win_sync(m_win);
    mpi::barrier_on_node();
    int disp_unit;
    MPI_Aint alloc_size;
    /*
     * MPI_Aint window_size0; int window_unit; double *win0_addr;
     * MPI_Win_shared_query(node_window, 0, &window_size0, &window_unit, &win0_addr);
     */
    ierr = MPI_Win_shared_query(m_win, 0, &alloc_size, &disp_unit, &m_data);
    REQUIRE_EQ(ierr, MPI_SUCCESS, "MPI Shared memory error");
    REQUIRE_EQ(uint_t(disp_unit), m_element_size, "incorrect window element size");
    REQUIRE_EQ(uint_t(alloc_size), m_nbyte, "incorrect total window size");
    MPI_Win_unlock_all(m_win);
    if (mpi::on_node_i_am_root()) std::memset(m_data, 0, m_nbyte);
    mpi::barrier_on_node();
}

SharedArrayBase::SharedArrayBase(SharedArrayBase&& rhs) :
        m_nelement(rhs.m_nelement), m_element_size(rhs.m_element_size), m_nbyte(rhs.m_nbyte){
    m_data = rhs.m_data;
    m_win = rhs.m_win;
    /*
     * nullify memory window handle in rhs so that the destructor does not
     * free memory still in use here
     */
    rhs.m_win = MPI_WIN_NULL;
}

SharedArrayBase::~SharedArrayBase() {
    REQUIRE_TRUE(m_data, "data pointer should be non-null");
    if (m_win != MPI_WIN_NULL) MPI_Win_free(&m_win);
}
