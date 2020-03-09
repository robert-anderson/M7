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
public:
    TableArray <table_T> m_send;
    table_T m_recv;
    template <typename factory_T>
    TableCommunicator(factory_T factory, size_t nrow) :
        m_send(factory, 4, nrow), m_recv(factory()) {
        m_recv.grow(nrow*mpi::nrank());
    }

    bool communicate() {
        defs::inds sendcounts(m_send.high_water_marks());
        for (auto &i : sendcounts) i *= m_send.row_length();
        defs::inds recvcounts(mpi::nrank(), 0ul);
        mpi::all_to_all(sendcounts, recvcounts);
        auto senddispls = m_send.offsets();
        defs::inds recvdispls(mpi::nrank(), 0ul);
        for (size_t i=1ul; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + sendcounts[i - 1];

        auto tmp = mpi::all_to_allv(m_send.data(), sendcounts, senddispls,
                                    m_recv.data(), recvcounts, recvdispls);

        m_recv.high_water_mark((recvdispls.back() + recvcounts.back()) / m_recv.row_length());
        m_send.zero();
        return tmp;
    }
};


#endif //M7_TABLECOMMUNICATOR_H
