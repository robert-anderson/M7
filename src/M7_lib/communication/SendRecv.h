//
// Created by rja on 01/10/22.
//

#ifndef M7_SENDRECV_H
#define M7_SENDRECV_H

#include "M7_lib/io/Logging.h"
#include "M7_lib/util/String.h"

#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/table/BufferedTableArray.h"

/**
 * A container for the send table array and recv table. Each element of the
 * send array corresponds to the "destination" rank of the MPI AllToAllV
 * communication invoked in the communicate method.
 * @tparam row_t
 *  Derived type of Row defining the data layout of both send and recv tables
 * @tparam mapped
 *  The send table is optionally mapped allowing rows with the same value in the
 *  mapped field to be accumulated together instead of occupying separate rows.
 *  This is a trade-off at the expense of more costly access (via hash tables)
 */
template<typename row_t, typename send_table_t>
class SendRecv {
    static_assert(std::is_base_of<Table<row_t>, send_table_t>::value, "Template args incompatible");
public:
    typedef BufferedTableArray<row_t, send_table_t> send_t;
    typedef buffered::Table<row_t> recv_t;
    typedef Table<row_t> recv_table_t;
private:
    send_t m_send;
    recv_t m_recv;

public:
    const row_t& m_row;
    /*
     * number of rows sent and recvd in the last call to communicate()
     * retained for stats reasons
     */
    uintv_t m_last_send_counts;
    uint_t m_last_recv_count = 0ul;

    /**
     * @param name
     *  base name for the buffers to be created for the send and recv tables
     * @param send
     *  send table instance, this is passed rather than a row, since this class allows for the definition of the send
     *  table as a MappedTable, whose ctor takes additional args
     * @param sizing
     *  buffer (re)sizing behaviour for both send and recv tables. estimated number of rows is taken to be in total
     *  across all nrank processes, and since there are nrank send tables, this estimate is divided by nrank to arrive
     *  at the estimated number of rows required in the send table
     */
    SendRecv(str_t name, const row_t &row, Sizing sizing) :
            m_send(name + " send", mpi::nrank(), row),
            m_recv(name + " recv", row), m_row(m_recv.m_row) {
        logging::info("Initially allocating {} per rank for each communicating buffer of \"{}\" (send and recv)",
                      string::memsize(mpi::nrank() * sizing.m_nrec_est * m_recv.row_size()), name);
        m_send.resize(integer::divceil(sizing.m_nrec_est, mpi::nrank()), 0.0);
        m_recv.resize(sizing.m_nrec_est, 0.0);
        m_send.set_expansion_factor(sizing.m_exp_fac);
        m_recv.set_expansion_factor(sizing.m_exp_fac);
    }

    uint_t row_size() const {
        return static_cast<const TableBase &>(m_recv).row_size();
    }

    send_t &send() {
        return m_send;
    }

    const send_t &send() const {
        return m_send;
    }

    typename send_t::table_t &send(const uint_t &i) {
        return m_send[i];
    }

    const typename send_t::table_t &send(const uint_t &i) const {
        return m_send[i];
    }

    Table<row_t> &recv() {
        return m_recv;
    }

    const Table<row_t> &recv() const {
        return m_recv;
    }

    void resize(uint_t nrow, double factor = -1.0) {
        m_send.resize(nrow, factor);
        m_recv.resize(nrow * mpi::nrank(), factor);
    }

    void expand(uint_t nrow, double factor = -1.0) {
        m_send.expand(nrow, factor);
        m_recv.expand(nrow * mpi::nrank(), factor);
    }

    /**
     * perform MPI alltoallv communication of the contents of all send buffers to all recv buffers then clear the send
     * table.
     */
    void communicate() {
        m_last_send_counts = m_send.hwms();
        uintv_t sendcounts(m_last_send_counts);
        for (auto &it: sendcounts) it *= row_size();
        uintv_t recvcounts(mpi::nrank(), 0ul);

//        auto all_sends_empty = !std::accumulate(sendcounts.cbegin(), sendcounts.cend(), 0ul);
//        if (all_sends_empty && m_last_recv_count) {
//            logging::debug_("this rank is sending no data at all, but it received data in the previous communication");
//        }
        m_recv.clear();

        mpi::all_to_all(sendcounts, recvcounts);

        auto senddispls = m_send.displs();
        uintv_t recvdispls(mpi::nrank(), 0ul);
        for (uint_t i = 1ul; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
        auto recv_size = recvdispls.back() + recvcounts.back();
        m_last_recv_count = recv_size / row_size();

        if (recv_size > static_cast<const TableBase &>(recv()).bw_size()) {
            /*
             * the recv table is full
             * this expansion by a factor is done explicitly here, because we
             * want to expand relative to size of the incoming data, not the
             * current size of the buffer
             */
            m_recv.resize(m_last_recv_count);
        }
        DEBUG_ASSERT_LE_ALL(recv_size, static_cast<const TableBase &>(recv()).bw_size(), "resize failed");

        REQUIRE_TRUE_ALL(m_send.begin(), "Send buffer is not allocated on all ranks!");
        REQUIRE_TRUE_ALL(m_recv.begin(), "Recv buffer is not allocated on all ranks!");

        auto tmp = mpi::all_to_allv(m_send.begin(), sendcounts, senddispls,
                                    m_recv.begin(), recvcounts, recvdispls);
        /*
         * check that the data addressed to this rank from this rank has been copied correctly
         */
        ASSERT(!send(mpi::irank()).begin() ||
               std::memcmp(send(mpi::irank()).begin(),
                           recv().begin() + recvdispls[mpi::irank()], recvcounts[mpi::irank()]) == 0);

        REQUIRE_TRUE_ALL(tmp, "MPI AllToAllV failed");
        recv().m_hwm = m_last_recv_count;
        m_send.clear();
    }

    void set_expansion_factor(double f) {
        m_send.set_expansion_factor(f);
        m_recv.set_expansion_factor(f);
    }
};

namespace send_recv {
    template <typename row_t> using BasicSend = SendRecv<row_t, Table<row_t>>;
    template <typename row_t> using MappedSend = SendRecv<row_t, MappedTable<row_t>>;
}


#endif //M7_SENDRECV_H
