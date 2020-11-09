//
// Created by rja on 09/11/2020.
//

#ifndef M7_COMMUNICATINGPAIR_H
#define M7_COMMUNICATINGPAIR_H

#include "BufferedTableArray.h"
#include "BufferedTable.h"
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/io/Logging.h"

template <typename table_t>
class CommunicatingPair {

    BufferedTableArray<table_t> m_send;
    BufferedTable<table_t> m_recv;

    template<typename ...Args>
    CommunicatingPair(Args... args): m_send(mpi::nrank(), args...), m_recv(args...){}

    size_t row_dsize() const {
        return static_cast<const TableX&>(m_recv).m_row_dsize;
    }

    void communicate() {
        auto hwms = m_send.hwms();
        defs::inds sendcounts(hwms);
        for (auto& it: sendcounts) it*=row_dsize();
        defs::inds recvcounts(mpi::nrank(), 0ul);

        //std::cout << "Sending elements " << utils::to_string(sendcounts) << std::endl;
        mpi::all_to_all(sendcounts, recvcounts);
        //std::cout << "Receiving elements " << utils::to_string(recvcounts) << std::endl;

        auto senddispls = m_send.displs();
        defs::inds recvdispls(mpi::nrank(), 0ul);
        for (size_t i = 1ul; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
       // std::cout << "Receiving displacements " << utils::to_string(recvdispls) << std::endl;

//    const auto ndword_send_tot = senddispls[mpi::nrank() - 1] + sendcounts[mpi::nrank() - 1];
//    const auto ndword_recv_tot = recvdispls[mpi::nrank() - 1] + recvcounts[mpi::nrank() - 1];


        logger::write("Send List usage fraction: " +
                      std::to_string(sendcounts[mpi::irank()] / double(m_send[0].buffer_dsize())), 0, logger::debug);
        logger::write("Receive List usage fraction: " +
                      std::to_string(recvcounts[mpi::irank()] / double(m_recv.buffer_dsize())), 0, logger::debug);

        ASSERT(recvcounts[mpi::irank()] < m_recv.dsize())
//    if (ndword_recv_tot > m_recv->dsize()) {
//        logger::write("Insufficient space for received data: reallocating recv list...");
//        m_recv->expand(ndword_recv_tot/m_recv->m_padded_row_dsize);
//    }

        auto tmp = mpi::all_to_allv(m_send.ptr(), sendcounts, senddispls,
                                    m_recv.ptr(), recvcounts, recvdispls);

        if (!tmp) throw std::runtime_error("MPI AllToAllV failed");

        static_cast<TableX&>(m_recv).m_hwm = (recvdispls.back() + recvcounts.back()) / row_dsize();
        std::cout << "Number of recvd elements " << m_recv->m_high_water_mark[0] << std::endl;
        m_send.clear();
        m_recv.clear();
    }
};


#endif //M7_COMMUNICATINGPAIR_H
