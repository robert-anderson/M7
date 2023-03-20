//
// Created by anderson on 20/03/2023.
//

#include "Reduction.h"

void reduction::Base::to_global_send() {
    auto& send_buf = mpi::g_send_reduce_buffers[m_itype];
    m_global_reduced_offset = send_buf.size();
    auto ptr = m_local_base.cbegin();
    send_buf.insert(send_buf.cend(), ptr, ptr + m_local_base.m_size);
    auto& recv_buf = mpi::g_recv_reduce_buffers[m_itype];
    if (recv_buf.size() < send_buf.size()) recv_buf.resize(send_buf.size());
}

void reduction::Base::from_global_recv() {
    DEBUG_ASSERT_NE(m_global_reduced_offset, ~0ul, "global offset should have been set in a prior to_global_send call");
    auto& buf = mpi::g_recv_reduce_buffers[m_itype];
    const auto ptr = m_reduced_base.begin();
    std::memcpy(ptr, buf.data()+m_global_reduced_offset, m_reduced_base.m_size);
    m_global_reduced_offset = ~0ul;
}

void reduction::all_reduce(mpi::Op op) {
    all_reduce_one_type<char>(op);
    all_reduce_one_type<short int>(op);
    all_reduce_one_type<int>(op);
    all_reduce_one_type<long int>(op);
    all_reduce_one_type<long long int>(op);
    all_reduce_one_type<unsigned char>(op);
    all_reduce_one_type<unsigned short int>(op);
    all_reduce_one_type<unsigned int>(op);
    all_reduce_one_type<unsigned long int>(op);
    all_reduce_one_type<unsigned long long int>(op);
    all_reduce_one_type<float>(op);
    all_reduce_one_type<double>(op);
    all_reduce_one_type<long double>(op);
    all_reduce_one_type<std::complex<float>>(op);
    all_reduce_one_type<std::complex<double>>(op);
    all_reduce_one_type<std::complex<long double>>(op);
    all_reduce_one_type<bool>(op);
}

void reduction::all_reduce(v_t<reduction::Base *> &members, mpi::Op op) {
    for (auto& member: members) member->to_global_send();
    all_reduce(op);
    for (auto& member: members) {
        member->from_global_recv();
        member->post_reduce();
    }
    for (auto& send: mpi::g_send_reduce_buffers) send.clear();
    for (auto& recv: mpi::g_recv_reduce_buffers) recv.assign(recv.size(), 0);
}
