//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/io/Logging.h>
#include "List.h"
#include "omp.h"

List::List(size_t nsegment) : Table(nsegment), m_high_water_mark(nsegment, 0ul) {}

void List::recv(List *list) {
    if (!compatible_with(*list)) throw std::runtime_error("Receiving List is not compatible");
    m_recv = list;
}


const defs::inds &List::high_water_mark() const {
    return m_high_water_mark;
}

const size_t &List::high_water_mark(const size_t isegment) const {
    return m_high_water_mark[isegment];
}

size_t List::push(const size_t &isegment) {
    ASSERT(isegment < m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment]++;
    if (tmp >= m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

size_t List::push(const size_t &isegment, const size_t &nrow) {
    ASSERT(isegment < m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment] += nrow;
    if (tmp > m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp-nrow;
}

size_t List::expand_push(const size_t &isegment, const size_t &nrow, double factor) {
    ASSERT(!omp_get_level()); // convenient but NOT threadsafe
    ASSERT(factor>=1.0);
    if (m_high_water_mark[isegment] + nrow > nrow_per_segment()){
        resize(factor*double(m_high_water_mark[isegment] + nrow));
    }
    return push(isegment, nrow);
}

void List::zero() {
    // TODO: no need to memset zero here, only included initially for clarity in debugging
    Table::zero();
    m_high_water_mark.assign(m_nsegment, 0ul);
}

void List::communicate() {
    ASSERT(m_recv);
    defs::inds sendcounts(m_high_water_mark);
    for (auto &i : sendcounts) i *= m_padded_row_dsize;
    defs::inds recvcounts(mpi::nrank(), 0ul);

    std::cout << "Sending elements " << utils::to_string(sendcounts) << std::endl;
    mpi::all_to_all(sendcounts, recvcounts);
    std::cout << "Receiving elements " << utils::to_string(recvcounts) << std::endl;

    auto senddispls = m_segment_doffsets;
    defs::inds recvdispls(mpi::nrank(), 0ul);
    for (size_t i = 1ul; i < mpi::nrank(); ++i)
        recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
    std::cout << "Receiving displacements " << utils::to_string(recvdispls) << std::endl;

    const auto ndword_send_tot = senddispls[mpi::nrank() - 1] + sendcounts[mpi::nrank() - 1];
    const auto ndword_recv_tot = recvdispls[mpi::nrank() - 1] + recvcounts[mpi::nrank() - 1];


    logger::write("Send List usage fraction: "+
                  std::to_string(ndword_send_tot/double(m_segment_dsize)), 0, logger::debug);
    logger::write("Receive List usage fraction: "+
    std::to_string(ndword_send_tot/double(m_segment_dsize)), 0, logger::debug);

    if (ndword_recv_tot > m_recv->dsize()) {
        logger::write("Insufficient space for received data: reallocating recv list...");
        m_recv->expand(ndword_recv_tot/m_recv->m_padded_row_dsize);
    }

    auto tmp = mpi::all_to_allv(m_data.data(), sendcounts, senddispls,
                                m_recv->m_data.data(), recvcounts, recvdispls);

    if (!tmp) throw std::runtime_error("MPI AllToAllV failed");

    m_recv->m_high_water_mark[0] = (recvdispls.back() + recvcounts.back()) / m_padded_row_dsize;
    std::cout << "Number of recvd elements " << m_recv->m_high_water_mark[0] << std::endl;
    zero();
}

void List::all_gather(List &local){
    ASSERT(compatible_with(local));
    defs::inds recvcounts(mpi::nrank());
    defs::inds displs(mpi::nrank(), 0ul);
    size_t nsend = local.high_water_mark(0)*m_padded_row_dsize;
    mpi::all_gather(&nsend, 1, recvcounts.data(), 1);
    for (size_t i=1ul; i<mpi::nrank(); ++i) displs[i] = displs[i-1]+recvcounts[i-1];
    size_t nrow = (displs.back()+recvcounts.back())/m_padded_row_dsize;
    resize(nrow);
    mpi::all_gatherv(local.m_data.data(), nsend, m_data.data(), recvcounts, displs);
    m_high_water_mark[0] = nrow;
}


void List::expand(size_t delta_nrow) {
    Table::expand(delta_nrow);
    if (m_recv) m_recv->expand(mpi::nrank()*delta_nrow);
}

std::string List::to_string() {
    return Table::to_string(m_high_water_mark);
}
