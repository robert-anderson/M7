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

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/conf/Conf.h>
#include <M7_lib/util/String.h>

#include "BufferedTable.h"
#include "BufferedTableArray.h"
#include "RowProtector.h"

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
class CommunicatingPair {
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
    CommunicatingPair(str_t name, const send_table_t &send, Sizing sizing) :
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
        if (mpi::nrank()==1ul) return;
        uint_t iblock, irank_dst;
        while (update(iblock, irank_dst)) m_moves.push_back({iblock, irank_dst});
    }
};



class Distribution {

    uintv_t m_block_iranks;

    /**
     * number of blocks stored on this MPI rank
     */
    uint_t m_nblock_local = 0ul;

public:
    uint_t nblock() const {
        return m_block_iranks.size();
    }

    uint_t nblock_() const {
        return m_nblock_local;
    }

    const uintv_t& block_iranks() const {
        return m_block_iranks;
    }

    Distribution(size_t nblock) {
        REQUIRE_GE(nblock, mpi::nrank(), "number of blocks may not be less than the number of ranks");
        m_block_iranks.reserve(nblock);
        /*
         * initialize distribution evenly
         */
        for (uint_t irank=0ul; irank < mpi::nrank(); ++irank){
            m_nblock_local = mpi::evenly_shared_count(nblock);
            for (uint_t iblock = 0ul; iblock < m_nblock_local; ++iblock) {
                m_block_iranks.push_back(irank);
            }
        }
        DEBUG_ASSERT_EQ(m_block_iranks.size(), nblock, "error in initial block allocation");
    }

    void update(const Redistributor& redist) {
        for (auto move: redist.m_moves) {
            auto& irank = m_block_iranks[move.m_iblock];
            if (mpi::i_am(irank)) --m_nblock_local;
            irank = move.m_dst_irank;
            if (mpi::i_am(irank)) ++m_nblock_local;
        }
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

class DistribDependentBase {
    virtual void on_redistribution() = 0;
};

struct CommunicatorBase {
    mutable std::set<DistribDependentBase*> m_dependents;
    ~CommunicatorBase() {
        REQUIRE_TRUE_ALL(m_dependents.empty(), "Communicator should never be outlived by its dependents");
    }
};

class DistribDependent : DistribDependentBase {
    const CommunicatorBase& m_comm_base;
public:
    DistribDependent(const CommunicatorBase& comm_base): m_comm_base(comm_base){
        m_comm_base.m_dependents.insert(this);
    }
    ~DistribDependent(){
        m_comm_base.m_dependents.erase(this);
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
struct Communicator : CommunicatorBase {
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
    typedef CommunicatingPair<comm_row_t, mapped_comm> comm_pair_t;
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
    typedef CommunicatingPair<store_row_t, false> redist_pair_t;
    redist_pair_t m_redist;
    /**
     * row protection information must also be shared
     */
    typedef SingleFieldRow<field::Number<uint_t>> prot_level_row_t;
    typedef CommunicatingPair<prot_level_row_t> prot_level_pair_t;
    prot_level_pair_t m_prot_level;
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
    Communicator(str_t name, const store_table_t &store, Sizing store_sizing,
                 const send_table_t &send, Sizing comm_sizing, uint_t nblock_per_rank):
            m_store(name + " store", store), m_comm(name, send, comm_sizing),
            m_name(name), m_dist(nblock_per_rank*mpi::nrank()),
            m_redist(name+" redistributor", store, comm_sizing),
            m_prot_level(name+" protection level", {{}}, comm_sizing),
            m_block_work_figures(m_dist.nblock()) {
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

    uint_t nrow_estimate(uint_t crude, double fac) {
        return std::max(1ul, crude) * fac;
    }

    Communicator(str_t name, const store_table_t &store, uint_t store_nrow_crude_est,
                 const send_table_t &send, uint_t comm_nrow_crude_est,
                 const conf::Buffers &buf_opts, const conf::LoadBalancing &ra_opts) :
        Communicator(name,
                store, {nrow_estimate(store_nrow_crude_est, buf_opts.m_store_fac_init), buf_opts.m_store_exp_fac},
                send, {nrow_estimate(comm_nrow_crude_est, buf_opts.m_comm_fac_init), buf_opts.m_comm_exp_fac},
                ra_opts.m_nblock_per_rank * mpi::nrank()) {}

    virtual ~Communicator() {}

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
            /*
             * sum so that all ranks have access to the same work figure information
             */
            auto local = m_block_work_figures;
            mpi::all_sum(local.data(), m_block_work_figures.data(), local.size());
        }
        /*
         * compute the best possible redistribution of blocks
         */
        Redistributor redist(m_dist.block_iranks(), m_block_work_figures, mpi::nrank());
        m_dist.update(redist);
        {
            /*
             * loop over all rows,
             *  write those which no longer belong here into the communicating send table,
             *  clear them from the store
             */
            auto& row = m_store.m_row;
            for (row.restart(); row.in_range(); row.step()) {
                if (row.key_field().is_zero()) continue;
                auto irank_owner = irank(row.key_field());
                if (!mpi::i_am(irank_owner)) {
                    m_redist.send(irank_owner).m_row.push_back_jump();
                    m_redist.send(irank_owner).m_row = row;
                    m_prot_level.send(irank_owner).m_row.push_back_jump();
                    m_prot_level.send(irank_owner).m_row.m_field = row.protection_level();
                    while (row.is_protected()) row.release();
                    m_store.clear_row(row);
                }
            }
            DEBUG_ASSERT_FALSE(m_redist.send(mpi::irank()).m_hwm, "rank should not be sending to itself");
        }
        m_redist.communicate();
        m_prot_level.communicate();
        {
            auto& recv_row = m_redist.recv().m_row;
            auto& prot_level_row = m_prot_level.recv().m_row;
            prot_level_row.restart();
            for (recv_row.restart(); recv_row.in_range(); recv_row.step()) {
                auto irank_owner = irank(recv_row.key_field());
                DEBUG_ASSERT_TRUE(mpi::i_am(irank_owner), "recv_row sent to wrong rank!");
                m_store.insert(recv_row.key_field());
                m_store.m_row = recv_row;
                for (uint_t ilevel=0ul; ilevel < uint_t(prot_level_row.m_field); ++ilevel) m_store.m_row.protect();
                prot_level_row.step();
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

namespace shared_rows {
    /**
     * A set of rows which is copied to all ranks, responds correctly to rank reallocation and is protected from erasure
     */
    template<typename store_row_t>
    struct Set : public DistribDependent {
        /**
         * the communicator whose m_store member stores the definitive row values
         */
        typedef MappedTable<store_row_t> store_t;
        const store_t& m_src_store;
        /**
         * set of row indices stored on this rank
         */
        std::set<uint_t> m_irows;
        /**
         * name to use for this row set in detailed logging
         */
        str_t m_name;
        /**
         * all rows gathered
         */
        BufferedTable<store_row_t, true> m_all;
        /**
         * send / recv tables for gathering into m_all
         */
        BufferedTable<store_row_t> m_gather_send;
        BufferedTable<store_row_t> m_gather_recv;

    public:
        template<typename comm_row_t, bool mapped_comm = false>
        Set(str_t name, const Communicator<store_row_t, comm_row_t, mapped_comm>& src_comm, uintv_t irows = {}) :
                DistribDependent(src_comm), m_src_store(src_comm.m_store),
                m_name(name),
                m_all(name+" all rows", {src_comm.m_store.m_row}),
                m_gather_send(name+" gather send", {src_comm.m_store.m_row}),
                m_gather_recv(name+" gather recv", {src_comm.m_store.m_row}) {
            for (auto irow: irows) add_(irow);
            full_update();
        }

        virtual ~Set() {}

        /**
         * bring all rows into m_gather_recv
         */
        void all_gatherv() {
            m_gather_send.clear();
            auto& send_row = m_gather_send.m_row;
            auto& src_row = m_src_store.m_row;
            for (auto irow: m_irows) {
                src_row.jump(irow);
                send_row.push_back_jump();
                send_row.copy_in(src_row);
            }
            static_cast<TableBase&>(m_gather_recv).all_gatherv(m_gather_send);
        }

        /**
         * direct buffer copy, no rows have changed locations
         */
        void update() {
            all_gatherv();
            m_all.m_bw = m_gather_recv.m_bw;
        }

        void row_index_update() {
            REQUIRE_EQ_ALL(m_all.m_hwm, m_irows.size(),
                           "this kind of refresh is only possible when all rows have already been gathered");
            m_irows = {};
            for (m_all.m_row.restart(); m_all.m_row.in_range(); m_all.m_row.step()) {
                //auto& row = m_src_store.m_row;
                const auto result = m_src_store[m_all.m_row.key_field()];
                if (result) m_irows.insert(*result);
            }
        }

        void full_update() {
            all_gatherv();
            m_all.clear();
            auto& recv_row = m_gather_recv.m_row;
            for (recv_row.restart(); recv_row.in_range(); recv_row.step()) m_all.insert(recv_row);
        }

        uint_t nrow_() const {
            return m_irows.size();
        }

//        uint_t nrow() const {
//            return mpi::all_sum(nrow_());
//        }

        void clear() {
            for (const auto& irow: m_irows) static_cast<const TableBase&>(m_src_store).release(irow);
            m_irows.clear();
        }

        void add_(uint_t irow) {
            static_cast<const TableBase&>(m_src_store).protect(irow);
            m_irows.insert(irow);
        }

        void add_(const store_row_t& row) {
            add_(row.index());
        }

    private:
        void on_redistribution() override {
            row_index_update();
            full_update();
        }
    };

    template<typename store_row_t>
    struct Single : Set<store_row_t> {
        using Set<store_row_t>::m_all;
        using Set<store_row_t>::clear;
        using Set<store_row_t>::add_;
        using Set<store_row_t>::full_update;

        template<typename comm_row_t, bool mapped_comm = false>
        Single(str_t name, const Communicator<store_row_t, comm_row_t, mapped_comm>& src_comm, TableBase::Loc loc):
            Set<store_row_t>(name, src_comm, loc) {
            m_all.m_row.restart();
        }

        void redefine(TableBase::Loc loc) {
            clear();
            if (loc.is_mine()) add_(loc.m_irow);
            full_update();
        }
    };
}


#endif //M7_COMMUNICATOR_H
