//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_COMMUNICATIONBUFFER_H
#define M7_COMMUNICATIONBUFFER_H

/*
class CommunicationBuffer : DataTable {

    void DataSystem::communicate() {
        MPIWrapper mpi;
        defs::inds sendcounts = m_send_buffer.highwatermark();
        for (auto &i : sendcounts) i*=m_send_buffer.total_datawords_used();
        defs::inds recvcounts(mpi.nrank(), 0ul);
        mpi.all_to_all(sendcounts, recvcounts);
        auto &senddispls = m_send_buffer.segment_dataword_offsets();
        //place the received data contiguously in the recv buffer:
        defs::inds recvdispls(mpi.nrank(), 0ul);
        for (auto i{1ul}; i<mpi.nrank(); ++i)
            recvdispls[i] = recvdispls[i-1]+sendcounts[i-1];

        mpi.all_to_allv(m_send_buffer.baseptr(), sendcounts, senddispls,
                        m_recv_buffer.baseptr(), recvcounts, recvdispls);
    }

};
*/

#endif //M7_COMMUNICATIONBUFFER_H
