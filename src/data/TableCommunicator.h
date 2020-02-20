//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_TABLECOMMUNICATOR_H
#define M7_TABLECOMMUNICATOR_H

#include <src/defs.h>
#include <src/parallel/MPIWrapper.h>
#include "TableArray.h"
#include "Specification.h"

template<typename table_T>
class TableCommunicator {
    TableArray <table_T> m_send;
    table_T m_recv;
public:
    TableCommunicator(Specification spec, size_t nrow_send, size_t nrow_recv) :
        m_send(mpi::nrank(), spec, nrow_send), m_recv(spec, nrow_recv) {}

    bool communicate(Table &recv) {
        defs::inds sendcounts(m_send.high_water_marks());
        for (auto &i : sendcounts) i *= m_send.row_length();
        defs::inds recvcounts(mpi::nrank(), 0ul);
        mpi::all_to_all(sendcounts, recvcounts);
        auto senddispls = m_send.offsets();
        defs::inds recvdispls(mpi::nrank(), 0ul);
        for (auto i{1ul}; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + sendcounts[i - 1];

        auto tmp = mpi::all_to_allv(m_send.data(), sendcounts, senddispls,
                                    recv.data(), recvcounts, recvdispls);
        //recv.set_highwatermark(0, (recvdispls.back() + recvcounts.back()) / row_length());
        //zero();
        return tmp;
    }

};


#endif //M7_TABLECOMMUNICATOR_H
