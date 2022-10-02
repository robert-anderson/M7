//
// Created by rja on 24/07/22.
//
//
// Created by Robert J. Anderson on 09/11/2020.
//

#ifndef M7_COMMUNICATOR_H
#define M7_COMMUNICATOR_H

#include <set>
#include <numeric>

#include "DistributedTable.h"

/**
 * Combines the SendRecv with a persistent storage table
 * @tparam store_row_t
 *  The Row class-derived data layout of the storage table
 * @tparam send_recv_row_t
 *  The Row class-derived data layout of the tables within the SendRecv object
 * @tparam mapped_comm
 *  optional mapping of send table in SendRecv
 */
template<typename store_row_t, typename send_recv_row_t, typename send_recv_table_t>
struct Communicator {
    static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
    static_assert(std::is_base_of<Row, send_recv_row_t>::value, "Template arg must be derived from Row");
    /**
     * the store table is a MappedTable which forms the persistent storage of the distributed object
     */
    typedef buffered::DistributedTable<store_row_t> store_t;
    store_t m_store;
    /**
     * convenient access to the distribution information (owned by the store table)
     */
    const Distribution& m_dist;
    /**
     * the communication is done by a pair of tables, send and recv. The Row type is not necessarily the same as that
     * used by the store table.
     */
    typedef SendRecv<send_recv_row_t, send_recv_table_t> send_recv_t;
    typedef typename send_recv_t::send_t::table_t send_table_t;
    send_recv_t m_send_recv;
    /**
     * identifier for the purposes of detailed logging (reallocations and remappings)
     */
    const str_t m_name;

    /**
     * @param name
     *  base name for all buffers created here
     * @param store
     *  store table instance
     * @param store_sizing
     *  buffer (re)sizing behaviour for the store table
     * @param send
     *  send table instance
     * @param comm_sizing
     *  buffer (re)sizing behaviour for the communicating pair
     */
    Communicator(str_t name,
                 const store_row_t &store_row, Sizing store_sizing,
                 const send_recv_row_t &send_recv_row, Sizing comm_sizing):
            m_store(store_row), m_dist(m_store.m_dist), m_send_recv(name, send, comm_sizing), m_name(name) {
        m_store.rename(name + " store");
        logging::info("Initially allocating {} per rank for store buffer of \"{}\" ",
                      string::memsize(store_sizing.m_nrow_est* m_store.record_size()), m_name);
        m_store.resize(store_sizing.m_nrow_est, 0.0);
        m_store.set_expansion_factor(store_sizing.m_exp_fac);
    }
    /**
     * ctor which uses configuration records
     * @param name
     *  base name for all buffers created here
     * @param store_nrow_crude_est
     *  "crude" estimate of the number of rows required in the store table. In general the max number of rows required
     *  cannot be known a-priori, so the best we can do is provide a sensible, order-of-magnitude estimation and allow
     *  the user to tune this via the initial scale factors in the buffer options sections of the input file
     * @param comm_nrow_crude_est
     *  "crude" estimate of the number of rows required per rank in the send table of the SendRecv
     * @param buf_opts
     *  options relating to the buffer allocation behavior for this object
     * @param ra_opts
     *  options relating to the RankAllocator behavior
     * @param store
     *  store table instance
     * @param send
     *  send table instance
     */

    uint_t nrow_estimate(uint_t crude, double fac) {
        return std::max(1ul, crude) * fac;
    }

//    Communicator(str_t name, const store_table_t &store, uint_t store_nrow_crude_est,
//                 const send_table_t &send, uint_t comm_nrow_crude_est,
//                 const conf::Buffers &buf_opts, const conf::LoadBalancing &ra_opts) :
//        Communicator(name,
//                store, {nrow_estimate(store_nrow_crude_est, buf_opts.m_store_fac_init), buf_opts.m_store_exp_fac},
//                send, {nrow_estimate(comm_nrow_crude_est, buf_opts.m_comm_fac_init), buf_opts.m_comm_exp_fac},
//                ra_opts.m_nblock_per_rank * mpi::nrank()) {}

    virtual ~Communicator() {}

    typename send_recv_t::send_t &send() {
        return m_send_recv.send();
    }

    const typename send_recv_t::send_t &send() const {
        return m_send_recv.send();
    }

    send_table_t &send(const uint_t &i) {
        return m_send_recv.send(i);
    }

    const send_table_t &send(const uint_t &i) const {
        return m_send_recv.send(i);
    }

    Table<send_recv_row_t> &recv() {
        return m_send_recv.recv();
    }

    const Table<send_recv_row_t> &recv() const {
        return m_send_recv.recv();
    }

    void communicate() {
        m_send_recv.communicate();
    }
};

namespace communicator {
    template <typename store_row_t, typename send_recv_row_t>
    using BasicSend = Communicator<store_row_t, send_recv_row_t, Table<send_recv_row_t>>;
    template <typename store_row_t, typename send_recv_row_t>
    using MappedSend = Communicator<store_row_t, send_recv_row_t, MappedTable<send_recv_row_t>>;
}


#endif //M7_COMMUNICATOR_H
