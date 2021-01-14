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
    double m_buffer_expansion_factor;

public:
    template<typename ...Args>
    CommunicatingPair(std::string name, double buffer_expansion_factor, Args... args):
    m_send(name+" send", mpi::nrank(), args...), m_recv(name+" recv", args...),
    m_buffer_expansion_factor(buffer_expansion_factor){}

    size_t row_dsize() const {
        return static_cast<const Table&>(m_recv).m_row_dsize;
    }

    BufferedTableArray<table_t>& send() {
        return m_send;
    }

    const BufferedTableArray<table_t>& send() const {
        return m_send;
    }

    table_t& send(const size_t& i) {
        return static_cast<table_t&>(m_send[i]);
    }

    const table_t& send(const size_t& i) const {
        return static_cast<const table_t&>(m_send[i]);
    }

    table_t& recv() {
        return static_cast<table_t&>(m_recv);
    }

    const table_t& recv() const {
        return static_cast<const table_t&>(m_recv);
    }

    void resize(size_t nrow) {
        m_send.resize(nrow);
        /*
         * recv buffer is dynamically resized during communication
         */
    }

    void expand(size_t nrow) {
        m_send.expand(nrow);
    }

    void expand(){
        m_send.expand();
    }

    void communicate() {

        auto hwms = m_send.hwms();
        defs::inds sendcounts(hwms);
        for (auto& it: sendcounts) it*=row_dsize();
        defs::inds recvcounts(mpi::nrank(), 0ul);

        //std::cout << "Sending datawords " << utils::to_string(sendcounts) << std::endl;
        mpi::all_to_all(sendcounts, recvcounts);
        //std::cout << "Receiving datawords " << utils::to_string(recvcounts) << std::endl;

        auto senddispls = m_send.displs();
        //std::cout << "Sending displacements " << utils::to_string(senddispls) << std::endl;
        defs::inds recvdispls(mpi::nrank(), 0ul);
        for (size_t i = 1ul; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
        //std::cout << "Receiving displacements " << utils::to_string(recvdispls) << std::endl;

//        logger::write("Send List usage fraction: " +
//                      std::to_string(sendcounts[mpi::irank()] / double(m_send[0].bw_dsize())), 0, logger::debug);
//        logger::write("Receive List usage fraction: " +
//                      std::to_string(recvcounts[mpi::irank()] / double(m_recv.bw_dsize())), 0, logger::debug);

        auto recv_dsize = recvdispls.back() + recvcounts.back();
        auto recv_nrow = recv_dsize / row_dsize();
        if (recv_dsize > m_recv.bw_dsize()){
            /*
             * the recv table is full
             */
            m_recv.resize(std::ceil((1.0+m_buffer_expansion_factor)*recv_nrow));
        }

        auto tmp = mpi::all_to_allv(m_send.dbegin(), sendcounts, senddispls,
                                    m_recv.dbegin(), recvcounts, recvdispls);

        /*
         * check that the data addressed to this rank from this rank has been copied correctly
         */
        ASSERT(mpi::nrank()> 0 or std::memcmp(
               (void*) (send(mpi::irank()).begin()),
               (void*) (recv().begin()+recvdispls[mpi::irank()]),
               recvcounts[mpi::irank()]*defs::nbyte_data)==0);

        if (!tmp) throw std::runtime_error("MPI AllToAllV failed");

        recv().m_hwm = recv_nrow;
//        std::cout << "Number of recvd datawords " << (recvdispls.back() + recvcounts.back()) << std::endl;
//        std::cout << "Number of recvd elements " << recv().m_hwm << std::endl;

        m_send.clear();
    }

    void set_expansion_factor(double f){
        m_send.set_expansion_factor(f);
        m_recv.set_expansion_factor(f);
    }
};


#endif //M7_COMMUNICATINGPAIR_H
