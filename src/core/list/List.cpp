//
// Created by Robert John Anderson on 2020-03-31.
//

#include "List.h"

List::List(size_t nsegment) : Table(nsegment), m_high_water_mark(nsegment, 0ul){}

void List::recv(List *list) {
    if(!compatible_with(*list)) throw std::runtime_error("Receiving List is not compatible");
    m_recv = list;
}


const defs::inds &List::high_water_mark() const {
    return m_high_water_mark;
}

size_t List::push(const size_t &isegment) {
    assert(isegment<m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment]++;
    if (tmp >= m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

size_t List::push(const size_t &isegment, const size_t &nrow) {
    assert(isegment<m_nsegment);
    size_t tmp;
#pragma omp atomic capture
    tmp = m_high_water_mark[isegment] += nrow;
    if (tmp >= m_nrow_per_segment) throw std::runtime_error("Reached capacity of List");
    return tmp;
}

void List::zero() {
    // TODO: no need to memset zero here, only included initially for clarity in debugging
    Table::zero();
    m_high_water_mark.assign(0ul, m_nsegment);
}

void List::communicate() {
    assert(m_recv);
    defs::inds sendcounts(m_high_water_mark);
    for (auto &i : sendcounts) i *= m_padded_row_dsize;
    defs::inds recvcounts(mpi::nrank(), 0ul);
    mpi::all_to_all(sendcounts, recvcounts);

    auto senddispls = m_segment_doffsets;
    defs::inds recvdispls(mpi::nrank(), 0ul);
    for (size_t i=1ul; i < mpi::nrank(); ++i)
        recvdispls[i] = recvdispls[i - 1] + sendcounts[i - 1];

    auto tmp = mpi::all_to_allv(m_data.data(), sendcounts, senddispls,
                                m_recv->m_data.data(), recvcounts, recvdispls);

    if(!tmp) throw std::runtime_error("MPI AllToAllV failed");

    m_recv->m_high_water_mark[0] = (recvdispls.back() + recvcounts.back()) / m_padded_row_dsize;
    zero();
}

void List::expand(size_t delta_nrow) {
    Table::expand(delta_nrow);
    if (m_recv) m_recv->expand(delta_nrow);
}
