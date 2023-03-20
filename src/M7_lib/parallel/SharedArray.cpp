//
// Created by rja on 24/01/23.
//

#include "SharedArray.h"

SharedArrayBase::SharedArrayBase(uint_t element_size) : m_element_size(element_size){}

void SharedArrayBase::alloc(uint_t nelement, uint_t element_size, MPI_Win *win, void **data) {
    const auto nbyte = nelement * element_size;
    const auto local_nbyte = mpi::on_node_i_am_root() ? nbyte : 0ul;
    auto ierr = MPI_Win_allocate_shared(local_nbyte, element_size, MPI_INFO_NULL, mpi::g_node_comm, data, win);
    REQUIRE_EQ(ierr, MPI_SUCCESS, "MPI Shared memory error");
    DEBUG_ASSERT_TRUE(*data, "data pointer not set");
    MPI_Win_lock_all(0, *win);
    MPI_Win_sync(*win);
    mpi::barrier_on_node();
    int disp_unit;
    MPI_Aint alloc_size;
    ierr = MPI_Win_shared_query(*win, 0, &alloc_size, &disp_unit, data);
    REQUIRE_EQ(ierr, MPI_SUCCESS, "MPI Shared memory error");
    REQUIRE_EQ(uint_t(disp_unit), element_size, "incorrect window element size");
    REQUIRE_EQ(uint_t(alloc_size), nbyte, "incorrect total window size");
    MPI_Win_unlock_all(*win);
    if (mpi::on_node_i_am_root()) std::memset(*data, 0, nbyte);
    mpi::barrier_on_node();
}

void SharedArrayBase::free(MPI_Win *win, void **data) {
    if (*data) {
        REQUIRE_FALSE(*win == MPI_WIN_NULL, "window should not be null if data is non-null");
        MPI_Win_free(win);
    }
    // otherwise, assume that the resource is still in use (has been moved)
}

void SharedArrayBase::alloc(uint_t nelement) {
    m_nelement = nelement;
    m_nbyte = nelement * m_element_size;
    alloc(nelement, m_element_size, &m_win, reinterpret_cast<void**>(&m_data));
    DEBUG_ASSERT_TRUE(m_data, "data pointer not set");
}

void SharedArrayBase::free() {
    auto data = reinterpret_cast<void*>(m_data);
    free(&m_win, &data);
}

SharedArrayBase::SharedArrayBase(uint_t nelement, uint_t element_size) : SharedArrayBase(element_size) {
    alloc(nelement);
}

SharedArrayBase &SharedArrayBase::operator=(const SharedArrayBase &other) {
    if (m_nbyte < other.m_nbyte) {
        free();
        alloc(other.m_nelement);
    }
    if (mpi::on_node_i_am_root()) std::memcpy(m_data, other.m_data, m_nbyte);
    return *this;
}

SharedArrayBase &SharedArrayBase::operator=(SharedArrayBase &&other) {
    REQUIRE_EQ(m_element_size, other.m_element_size, "incompatible shared arrays");
    free();
    m_nelement = other.m_nelement;
    m_nbyte = other.m_nbyte;
    m_data = other.m_data;
    m_win = other.m_win;
    // nullify pointer in other so that the destructor does not free memory still in use here
    other.m_data = nullptr;
    return *this;
}

SharedArrayBase::SharedArrayBase(const SharedArrayBase &other) : SharedArrayBase(other.m_nelement, other.m_element_size) {
    *this = other;
}

SharedArrayBase::SharedArrayBase(SharedArrayBase &&other) : SharedArrayBase(other.m_element_size) {
    *this = std::move(other);
}

SharedArrayBase::~SharedArrayBase() {
    free();
}