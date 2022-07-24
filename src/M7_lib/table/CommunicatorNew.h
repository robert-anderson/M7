//
// Created by rja on 24/07/22.
//
//
// Created by Robert J. Anderson on 09/11/2020.
//

#ifndef M7_COMMUNICATOR_NEW_H
#define M7_COMMUNICATOR_NEW_H

#include <set>
#include <numeric>

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/conf/Conf.h>
#include <M7_lib/util/String.h>

#include "BufferedTable.h"
#include "BufferedTableArray.h"
#include "RowProtector.h"
#include "M7_lib/parallel/Redistributor.h"

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
template<typename row_t, bool mapped = false>
class CommunicatingPairNew {
public:
    typedef BufferedTableArray<row_t, mapped> send_t;
    typedef BufferedTable<row_t> recv_t;
    typedef typename send_t::table_t send_table_t;
    typedef typename recv_t::table_t recv_table_t;
private:
    send_t m_send;
    recv_t m_recv;

public:
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
     *  buffer (re)sizing behaviour for both send and recv tables
     */
    CommunicatingPairNew(str_t name, const send_table_t &send, Sizing sizing) :
            m_send(name + " send", mpi::nrank(), send),
            m_recv(name + " recv", recv_table_t(send.m_row)) {
        logging::info("Initially allocating {} per rank for each communicating buffer of \"{}\" (send and recv)",
                      string::memsize(mpi::nrank() * sizing.m_nrow_est*m_recv.row_size()), name);
        m_send.resize(sizing.m_nrow_est, 0.0);
        m_recv.resize(mpi::nrank() * sizing.m_nrow_est, 0.0);
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
        m_send.expand(nrow);
        m_recv.expand(nrow * mpi::nrank());
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

        auto all_sends_empty = !std::accumulate(sendcounts.cbegin(), sendcounts.cend(), 0ul);
        DEBUG_ONLY(all_sends_empty);
        if (all_sends_empty && m_last_recv_count) {
            logging::debug_("this rank is sending no data at all, but it received data in the previous communication");
        }

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
        ASSERT(!send(mpi::irank()).begin() or
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

/**
 * decide which block transfers to make in order to equalize load across all ranks
 */
class Redistributor {
    /**
     * map from block index to rank index
     */
    const uintv_t& m_block_iranks;
    /**
     * workload figure for each block
     */
    const v_t<double>& m_block_work_figs;
    /**
     * map from rank index to block indices
     */
    v_t<std::set<uint_t>> m_rank_iblocks;
    /**
     * total workload figure for each rank
     */
    v_t<double> m_total_work_figs;
    /**
     * total workload figure per rank with perfect load balancing
     */
    const double m_perfect_work_fig;
    /**
     * flag set to true if the rank is busier than or equal to the average load (for checks)
     */
    const v_t<bool> m_busier_or_avg;

    v_t<bool> make_busier_or_avg() const {
        v_t<bool> flags;
        flags.reserve(m_total_work_figs.size());
        for (auto v: m_total_work_figs) flags.push_back(v>=m_perfect_work_fig);
        return flags;
    }

    /**
     * @param irank_busiest
     *  index of busiest rank (after already decided transfers)
     * @param limit
     *  maximum value of workload figure
     * @return
     *  index of block with the largest workload figure less than limit
     */
    uint_t iblock_max_lt(uint_t irank_busiest, double limit) {
        uint_t iblock_best = ~0ul;
        for (auto iblock: m_rank_iblocks[irank_busiest]) {
            if ((m_block_work_figs[iblock] < limit) && (iblock_best==~0ul ||
                (m_block_work_figs[iblock] > m_block_work_figs[iblock_best]))) iblock_best = iblock;
        }
        return iblock_best;
    }

    /**
     * do one iteration of the reallocation algorithm
     * @param iblock
     *  index of the block to transfer
     * @param irank_dst
     *  index of MPI rank to which block iblock is destined
     * @return
     *  false if no further transfers could be found
     */
    bool update(uint_t& iblock, uint_t& irank_dst) {
        auto it_src = std::max_element(m_total_work_figs.begin(), m_total_work_figs.end());
        const uint_t irank_src = std::distance(m_total_work_figs.begin(), it_src);
        auto it_dst = std::min_element(m_total_work_figs.begin(), m_total_work_figs.end());
        irank_dst = std::distance(m_total_work_figs.begin(), it_dst);
        DEBUG_ASSERT_NE(irank_src, irank_dst, "rank should never transfer to itself");

        // edge case where one of the extreme ranks has exactly the perfect load
        if (*it_src == m_perfect_work_fig || *it_dst == m_perfect_work_fig) {
            DEBUG_ASSERT_TRUE(*it_src == m_perfect_work_fig || *it_dst == m_perfect_work_fig,
                              "both should be perfectly distributed");
            return false;
        }

        DEBUG_ASSERT_GT(*it_src, m_perfect_work_fig,
                        "lazier than average rank should never be identified as busiest");
        DEBUG_ASSERT_LT(*it_dst, m_perfect_work_fig,
                        "busier than average rank should never be identified as laziest");

        /*
         * find the index of the block with the largest weight less than half the difference between the work figures of
         * the current busiest and laziest ranks. this will be the most efficient move. if the work figure is greater
         * than this limit, the lazy rank will become busier than average
         */
        iblock = iblock_max_lt(irank_src, (*it_src - *it_dst) / 2.0);
        if (iblock<m_block_iranks.size()) {
            DEBUG_ASSERT_FALSE(m_rank_iblocks[irank_src].empty(), "source rank should have at least one block");
            const auto diff = m_block_work_figs[iblock];
            *it_src -= diff;
            *it_dst += diff;
            m_rank_iblocks[irank_src].erase(iblock);
            m_rank_iblocks[irank_dst].insert(iblock);
            return true;
        }
        return false;
    }

public:

    struct Move {
        uint_t m_iblock, m_dst_irank;
        bool operator==(const Move& other) const {
            return m_iblock==other.m_iblock && m_dst_irank==other.m_dst_irank;
        }
    };
    v_t<Move> m_moves;

    Redistributor(const uintv_t& block_iranks, const v_t<double>& block_work_figs, uint_t nrank):
        m_block_iranks(block_iranks), m_block_work_figs(block_work_figs), m_rank_iblocks(nrank),
        m_perfect_work_fig(std::accumulate(m_block_work_figs.cbegin(), m_block_work_figs.cend(), 0.0)/nrank),
        m_busier_or_avg(make_busier_or_avg()){

        {
            /*
             * construct the inverse one-to-many map from ranks to blocks
             */
            uint_t iblock = 0ul;
            for (auto& irank: block_iranks) m_rank_iblocks[irank].insert(iblock++);
        }

        m_total_work_figs.resize(nrank);
        for (uint_t iblock=0ul; iblock<m_block_iranks.size(); ++iblock) {
            m_total_work_figs[m_block_iranks[iblock]] += m_block_work_figs[iblock];
            DEBUG_ASSERT_GE(m_block_work_figs[iblock], 0.0, "block work figure should be non-negative");
        }
        uint_t iblock, irank_dst;
        while (update(iblock, irank_dst)) m_moves.push_back({iblock, irank_dst});
    }
};



class Distribution {

    const uint_t m_nrank;
    uintv_t m_block_iranks;

public:
    uint_t nblock() const {
        return m_block_iranks.size();
    }

    const uintv_t& block_iranks() const {
        return m_block_iranks;
    }

    Distribution(size_t nblock, size_t nrank): m_nrank(nrank) {
        REQUIRE_GE(nblock, nrank, "number of blocks may not be less than the number of ranks");
        m_block_iranks.reserve(nblock);
        /*
         * initialize distribution evenly
         */
        for (uint_t irank=0ul; irank < nrank; ++irank){
            for (uint_t iblock = 0ul; iblock < uint_t(mpi::evenly_shared_count(nblock)); ++iblock)
                m_block_iranks.push_back(irank);
        }
        DEBUG_ASSERT_EQ(m_block_iranks.size(), nblock, "error in initial block allocation");
    }

    void update(const Redistributor& redist) {
        for (auto move: redist.m_moves) m_block_iranks[move.m_iblock] = move.m_dst_irank;
    }

    template<typename field_t>
    uint_t iblock(const field_t& field) const {
        return field.hash()%m_block_iranks.size();
    }

    template<typename field_t>
    uint_t irank(const field_t& field) const {
        return m_block_iranks[iblock(field)];
    }
};




/**
 * Combines the CommunicatingPair with a persistent storage table
 * @tparam store_row_t
 *  The Row class-derived data layout of the storage table
 * @tparam comm_row_t
 *  The Row class-derived data layout of the tables within the CommunicatingPair
 * @tparam mapped_comm
 *  optional mapping of send table in CommunicatingPair
 */
template<typename store_row_t, typename comm_row_t, bool mapped_comm = false>
struct CommunicatorNew {
    static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
    static_assert(std::is_base_of<Row, comm_row_t>::value, "Template arg must be derived from Row");
    /**
     * the store table is a MappedTable which forms the persistent storage of the distributed object
     */
    typedef BufferedTable<store_row_t, true> store_t;
    typedef typename store_t::table_t store_table_t;
    store_t m_store;
    /**
     * the key field of the store table determines the block, and therefore MPI rank index, to which the row belongs
     */
    typedef typename KeyField<store_row_t>::type key_field_t;
    /**
     * the communication is done by a pair of tables, send and recv. The Row type is not necessarily the same as that
     * used by the store table.
     */
    typedef CommunicatingPairNew<comm_row_t, mapped_comm> comm_pair_t;
    typedef typename comm_pair_t::send_t::table_t send_table_t;
    comm_pair_t m_comm;
    /**
     * identifier for the purposes of detailed logging (reallocations and remappings)
     */
    str_t m_name;
    /**
     * current allocation of load balancing blocks to MPI rank indices
     */
    Distribution m_dist;
    /**
     * redistribution of rows is crucial to achieve load balancing. this pair must always have the same Row type as the
     * store table though
     */
    typedef CommunicatingPairNew<store_row_t, false> redist_pair_t;
    redist_pair_t m_redist;
    /**
     * redistribution requires information about the amount of work done (by an arbitrary measure) by each block.
     */
    v_t<double> m_block_work_figures;

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
     * @param nblock_per_rank
     *  number of distribution blocks per rank
     */
    CommunicatorNew(str_t name, const store_table_t &store, Sizing store_sizing,
                 const send_table_t &send, Sizing comm_sizing, uint_t nblock_per_rank):
            m_store(name + " store", store), m_comm(name, send, comm_sizing),
            m_name(name), m_dist(nblock_per_rank*mpi::nrank(), mpi::nrank()),
            m_redist(name+" redistributor", store, comm_sizing), m_block_work_figures(m_dist.nblock()) {
        logging::info("Initially allocating {} per rank for store buffer of \"{}\" ",
                      string::memsize(store_sizing.m_nrow_est*m_store.row_size()), name);
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
     *  "crude" estimate of the number of rows required per rank in the send table of the CommunicatingPair
     * @param buf_opts
     *  options relating to the buffer allocation behavior for this object
     * @param ra_opts
     *  options relating to the RankAllocator behavior
     * @param store
     *  store table instance
     * @param send
     *  send table instance
     */
//    Communicator(str_t name, uint_t store_nrow_crude_est, uint_t comm_nrow_crude_est,
//                 const conf::Buffers &buf_opts, const conf::LoadBalancing &ra_opts,
//                 const store_table_t &store, const send_table_t &send) :
//            Communicator(name, std::max(1ul, store_nrow_crude_est) * buf_opts.m_store_fac_init, buf_opts.m_store_exp_fac,
//                         std::max(1ul, comm_nrow_crude_est) * buf_opts.m_comm_fac_init, buf_opts.m_comm_exp_fac,
//                         store, send,
//                         ra_opts.m_nblock_per_rank * mpi::nrank(), ra_opts.m_period, ra_opts.m_acceptable_imbalance,
//                         ra_opts.m_nnull_updates_deactivate) {}

    virtual ~CommunicatorNew() {}

    /**
     * @param key
     *  m_store-mapping key instance
     * @return
     *  rank index assigned to the key specified
     */
    uint_t irank(const key_field_t& key) const {
        return m_dist.irank(key);
    }

    void accumulate_work_figure(const key_field_t& key, double work_done) {
        m_block_work_figures[m_dist.iblock(key)]+=work_done;
    }

    void clear_work_figures() {
        m_block_work_figures.assign(m_block_work_figures.size(), 0.0);
    }

    void redistribute() {
        {
            auto local = m_block_work_figures;
            mpi::all_sum(local.data(), m_block_work_figures.data(), local.size());
        }
        Redistributor redist(m_dist.block_iranks(), m_block_work_figures, mpi::nrank());
        m_dist.update(redist);
        {
            auto& row = m_store.m_row;
            for (row.restart(); row.in_range(); row.step()) {
                auto irank_owner = irank(row.key_field());
                if (!mpi::i_am(irank_owner)) {
                    m_redist.send(irank_owner).m_row.push_back_jump();
                    m_redist.send(irank_owner).m_row = row;
                }
            }
            DEBUG_ASSERT_FALSE(m_redist.send(mpi::irank()).m_hwm, "rank should not be sending to itself");
        }
        m_redist.communicate();
        {
            auto& row = m_redist.recv().m_row;
            for (row.restart(); row.in_range(); row.step()) {
                auto irank_owner = irank(row.key_field());
                DEBUG_ASSERT_TRUE(mpi::i_am(irank_owner), "row sent to wrong rank!");
                m_store.insert(row.key_field());
                m_store.m_row = row;
            }
        }
        clear_work_figures();
    }

    typename comm_pair_t::send_t &send() {
        return m_comm.send();
    }

    const typename comm_pair_t::send_t &send() const {
        return m_comm.send();
    }

    send_table_t &send(const uint_t &i) {
        return m_comm.send(i);
    }

    const send_table_t &send(const uint_t &i) const {
        return m_comm.send(i);
    }

    Table<comm_row_t> &recv() {
        return m_comm.recv();
    }

    const Table<comm_row_t> &recv() const {
        return m_comm.recv();
    }

    void communicate() {
        m_comm.communicate();
    }
};

/**
 * Extends the Communicator with the introduction of redistribution
 * @tparam store_row_t
 *  The Row class-derived data layout of the storage table
 * @tparam comm_row_t
 *  The Row class-derived data layout of the tables within the CommunicatingPair
 * @tparam mapped_comm
 *  optional mapping of send table in CommunicatingPair
 */
template<typename store_row_t, typename comm_row_t, bool mapped_comm = false>
struct DynamicCommunicator : CommunicatorNew<store_row_t, comm_row_t, mapped_comm> {
    typedef CommunicatorNew<store_row_t, comm_row_t, mapped_comm> base_t;
    typedef typename base_t::store_table_t store_table_t;
    typedef typename base_t::send_table_t send_table_t;
    typedef typename base_t::key_field_t key_field_t;

    typedef CommunicatingPairNew<store_row_t, false> redist_pair_t;
    redist_pair_t m_redist;

    DynamicCommunicator(str_t name, const store_table_t &store, Sizing store_sizing,
                        const send_table_t &send, Sizing comm_sizing, uint_t nblock_per_rank):
            base_t(name, store, store_sizing, send, comm_sizing, nblock_per_rank){}

};


#endif //M7_COMMUNICATOR_NEW_H
