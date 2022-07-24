//
// Created by rja on 17/07/22.
//

#ifndef M7_REDISTRIBUTOR_H
#define M7_REDISTRIBUTOR_H

#include <numeric>
#include <algorithm>
#include "MPIAssert.h"
#include "M7_lib/table/BufferedTableArray.h"
#include "M7_lib/table/BufferedTable.h"
#include "RankAllocator.h"



struct LoadBalancerBase {
    /**
     * workload figures for each block
     */
    v_t<double> m_block_work_figs;

    uintv_t m_block_to_rank;

    void accumulate_work_fig(uint_t iblock, double cost){
        m_block_work_figs[iblock]+=cost;
    }
};

/*
 * template<typename store_row_t, typename comm_row_t, bool mapped_send = false>
struct Communicator {
    static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
    static_assert(std::is_base_of<Row, comm_row_t>::value, "Template arg must be derived from Row");

    typedef typename KeyField<store_row_t>::type key_field_t;
    typedef BufferedTable<store_row_t, true> store_t;
struct LoadBalancer {
    m_lb(name, m_store, nblock_ra, period_ra, acceptable_imbalance, nnull_updates_deactivate),

};
 */


/*
 * a load balancing exchange involves an all-to-all communication
 *
 */


namespace lb {

    /**
     * A container for the send table array and recv table. Each element of the
     * send array corresponds to the "destination" rank of the MPI AllToAllV
     * communication invoked in the communicate method.
     * @tparam row_t
     *  Derived type of Row defining the data layout of both send and recv tables
     * @tparam mapped_send
     *  The send table is optionally mapped allowing rows with the same value in the
     *  mapped field to be accumulated together instead of occupying separate rows.
     *  This is a trade-off at the expense of more costly access (via hash tables)
     */
    template<typename row_t, bool mapped_send = false>
    class CommunicatingPair {
    public:
        typedef BufferedTableArray<row_t, mapped_send> send_t;
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
         *  buffer sizing/resizing options
         */
        CommunicatingPair(str_t name, const send_table_t &send, Sizing sizing) :
                m_send(name + " send", mpi::nrank(), send),
                m_recv(name + " recv", recv_table_t(send.m_row)) {
            logging::info("Initially allocating {} per rank for each communicating buffer of \"{}\" (send and recv)",
                          string::memsize(mpi::nrank() * sizing.m_nrow_est * m_recv.row_size()), name);
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
                logging::debug_(
                        "this rank is sending no data at all, but it received data in the previous communication");
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
     * Combines the CommunicatingPair with a persistent storage table
     * @tparam store_row_t
     *  The Row class-derived data layout of the storage table
     * @tparam comm_row_t
     *  The Row class-derived data layout of the tables within the CommunicatingPair
     * @tparam mapped_send
     *  optional mapping of sending tables of CommunicatingPair
     */
    template<typename store_row_t, typename comm_row_t, bool mapped_send = false>
    struct Communicator {
        static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
        static_assert(std::is_base_of<Row, comm_row_t>::value, "Template arg must be derived from Row");

        typedef typename KeyField<store_row_t>::type key_field_t;
        typedef BufferedTable<store_row_t, true> store_t;
        typedef CommunicatingPair<comm_row_t, mapped_send> comm_t;
        typedef typename store_t::table_t store_table_t;
        typedef typename comm_t::send_t::table_t send_table_t;

        store_t m_store;
        comm_t m_comm;
        str_t m_name;

        /**
         * @param name
         *  base name for all buffers created here
         * @param store
         *  store table instance
         * @param store_sizing
         *  buffer sizing options for the store table
         * @param send
         *  send table instance defining the communicating pair (send/recv)
         * @param comm_sizing
         *  buffer sizing options for the communicating pair table
         */
        Communicator(str_t name, const store_table_t &store, Sizing store_sizing,
                     const send_table_t &send, Sizing comm_sizing) :
                m_store(name + " store", store), m_comm(name, send, comm_sizing), m_name(name) {
            logging::info("Initially allocating {} per rank for store buffer of \"{}\" ",
                          string::memsize(store_sizing.m_nrow_est * m_store.row_size()), name);
            m_store.resize(store_sizing.m_nrow_est, 0.0);
            m_store.set_expansion_factor(store_sizing.m_exp_fac);
        }
#if 0
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
        Communicator(str_t name, uint_t store_nrow_crude_est, uint_t comm_nrow_crude_est,
                     const conf::Buffers &buf_opts, const conf::LoadBalancing &ra_opts,
                     const store_table_t &store, const send_table_t &send) :
                Communicator(name, std::max(1ul, store_nrow_crude_est) * buf_opts.m_store_fac_init, buf_opts.m_store_exp_fac,
                             std::max(1ul, comm_nrow_crude_est) * buf_opts.m_comm_fac_init, buf_opts.m_comm_exp_fac,
                             store, send,
                             ra_opts.m_nblock_per_rank * mpi::nrank(), ra_opts.m_period, ra_opts.m_acceptable_imbalance,
                             ra_opts.m_nnull_updates_deactivate) {}
#endif

        virtual ~Communicator() {}

        typename comm_t::send_t &send() {
            return m_comm.send();
        }

        const typename comm_t::send_t &send() const {
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


#if 0
    /**
     * Combines the CommunicatingPair with a persistent storage table and a RankAllocator
     * @tparam store_row_t
     *  The Row class-derived data layout of the storage table
     * @tparam comm_row_t
     *  The Row class-derived data layout of the tables within the CommunicatingPair
     * @tparam mapped_send
     *  optional mapping of CommunicatingPair
     */
    template<typename store_row_t, typename comm_row_t, bool mapped_send = false>
    struct Communicator {
        static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
        static_assert(std::is_base_of<Row, comm_row_t>::value, "Template arg must be derived from Row");

        typedef typename KeyField<store_row_t>::type key_field_t;
        typedef BufferedTable<store_row_t, true> store_t;
        typedef CommunicatingPair<comm_row_t, mapped_send> comm_pair_t;
        typedef typename store_t::table_t store_table_t;
        typedef typename comm_pair_t::send_t::table_t send_table_t;

        store_t m_store;
        comm_pair_t m_comm;
        LoadBalancer m_lb;
        str_t m_name;

        /**
         * A set of rows which responds correctly to rank reallocation and is protected from erasure
         */
        struct DynamicRowSet : RankDynamic, RowProtector {
            typedef RankAllocator<store_row_t> ra_t;
            /**
             * the communicator whose m_store member stores the definitive row values
             */
            const Communicator &m_source;
            /**
             * set of dynamic row indices stored on this rank
             */
            std::set<uint_t> m_irows;
            /**
             * the counts and displs determine a global index for a local dynamic row. The ith element of m_irows on MPI
             * rank j has global index m_displs[j]+i
             */
            uintv_t m_counts;
            uintv_t m_displs;
            uintv_t m_ranks_with_any_rows;

            str_t m_name;
            /*
             * row = row in mapped table
             * trow = row among transferred rows
             * drow = dynamic row index among transferred rows
             */
            const int m_ntrow_to_track_p2p_tag;
            const int m_itrows_to_track_p2p_tag;
            uintv_t m_idrows;
            uint_t m_itrow;
            uint_t m_ndrow_found;

        public:
            DynamicRowSet(const Communicator &comm, str_t name) :
                    RankDynamic(comm.m_ra),
                    RowProtector(comm.m_store),
                    m_source(comm),
                    m_counts(mpi::nrank(), 0ul),
                    m_displs(mpi::nrank(), 0ul),
                    m_name(name),
                    m_ntrow_to_track_p2p_tag(mpi::new_p2p_tag()),
                    m_itrows_to_track_p2p_tag(mpi::new_p2p_tag()) {
                logging::debug("P2P tag for number of dynamic row indices to transfer for \"{}\": {}", m_name,
                               m_ntrow_to_track_p2p_tag);
                logging::debug("P2P tag for array of dynamic row indices to transfer for \"{}\": {}", m_name,
                               m_itrows_to_track_p2p_tag);
                m_ranks_with_any_rows.reserve(mpi::nrank());
            }

            virtual ~DynamicRowSet() {}

            uint_t nrow_() const {
                return m_irows.size();
            }

            uint_t nrow() const {
                return mpi::all_sum(nrow_());
            }

            void clear() {
                for (const auto &irow: m_irows) RowProtector::release(irow);
                m_irows.clear();
            }

            void add_(uint_t irow) {
                RowProtector::protect(irow);
                m_irows.insert(irow);
            }

            void add_(const store_row_t &row) {
                add_(row.index());
            }

        public:
            virtual void update() {
                mpi::all_gather(nrow_(), m_counts);
                mpi::counts_to_displs_consec(m_counts, m_displs);
                m_ranks_with_any_rows.clear();
                for (uint_t irank = 0; irank < mpi::nrank(); ++irank) {
                    if (m_counts[irank]) m_ranks_with_any_rows.push_back(irank);
                }
            }

            bool has_row(uint_t irow) override {
                return m_irows.find(irow) != m_irows.end();
            }

            void before_block_transfer(const uintv_t &irows_send, uint_t irank_send, uint_t irank_recv) override {
                uint_t nrow_transfer;
                if (mpi::i_am(irank_send)) {
                    m_idrows.clear();
                    // look for Dynamic rows among those being transferred
                    uint_t itrow = 0ul;
                    for (const auto &irow: irows_send) {
                        if (m_irows.find(irow) != m_irows.end()) {
                            /*
                             * this row is tracked by this DynamicRowSet, so the recving rank needs to know
                             * which of the incoming transferred rows it is now responsible for tracking!
                             */
                            m_idrows.push_back(itrow);
                            if (is_protected(irow)) release(irow);
                            m_irows.erase(irow);
                        }
                        ++itrow;
                    }
                    nrow_transfer = m_idrows.size();
                    mpi::send(&nrow_transfer, 1, irank_recv, m_ntrow_to_track_p2p_tag);
                    // only bother with vector send if there's anything to transfer
                    if (!nrow_transfer) {
                        logging::debug_("No dynamic rows to transfer for \"{}\"", m_name);
                        return;
                    }
                    logging::debug_("Sending {} dynamic rows for \"{}\" from rank {} to {}",
                                    nrow_transfer, m_name, irank_send, irank_recv);
                    mpi::send(m_idrows.data(), nrow_transfer, irank_recv, m_itrows_to_track_p2p_tag);
                }

                if (mpi::i_am(irank_recv)) {
                    // rezero this counter, it will be incremented in every call to on_row_recv_
                    m_itrow = 0ul;
                    // rezero this counter, it will be incremented whenever on_row_recv_ finds a dynamic row
                    m_ndrow_found = 0ul;

                    m_idrows.clear();
                    mpi::recv(&nrow_transfer, 1, irank_send, m_ntrow_to_track_p2p_tag);
                    if (!nrow_transfer) {
                        logging::debug_("No dynamic rows to transfer for \"{}\"", m_name);
                        return;
                    }
                    logging::debug_("Recving {} dynamic rows for \"{}\" to rank {} from {}",
                                    nrow_transfer, m_name, irank_recv, irank_send);
                    m_idrows.resize(nrow_transfer, ~0ul);
                    mpi::recv(m_idrows.data(), nrow_transfer, irank_send, m_itrows_to_track_p2p_tag);
                    logging::debug_("idrows: {}", convert::to_string(m_idrows));
                }
            }

            void on_row_recv_(uint_t irow) override {
                // check if there's any more rows we may need to track
                if (m_ndrow_found < m_idrows.size()) {
                    ASSERT(m_ndrow_found < m_idrows.size());
                    const auto next_itrow = m_idrows[m_ndrow_found];
                    if (m_itrow == next_itrow) {
                        // this transferred row is dynamic
                        DEBUG_ASSERT_TRUE(m_irows.find(irow) == m_irows.end(),
                                          "Transferred dynamic row shouldn't already be here!");
                        m_irows.insert(irow);
                        protect(irow);
                        ++m_ndrow_found;
                    }
                }
                ++m_itrow;
            }

            void after_block_transfer() override {
                update();
            }
        };

        /**
         * A set of dynamic rows from which a subset of fields can be loaded into a contiguous local table and communicated to a
         * contiguous global table. That is to say: only partially shared
         */
        template<typename contig_row_t>
        struct PartSharedRowSet : public DynamicRowSet {
            using DynamicRowSet::m_source;
            using DynamicRowSet::nrow;
            using DynamicRowSet::m_displs;
            using DynamicRowSet::m_counts;
            using DynamicRowSet::m_irows;
            using DynamicRowSet::m_name;
            /**
             * the (mapped) table which loads data from rows of m_source.m_store (arbitrary order, non-contiguous) into m_global
             * (contiguous).
             */
            BufferedTable<contig_row_t> m_local;
            /**
             * the table which holds data extracted from all tracked rows in the dynamic set. these rows copies are refreshed
             * with a call to update method. This table is intended for reading data from all MPI ranks, as such there is no
             * mechanism for committing changes to m_source.m_store from m_global. It is to be treated as a read-only view of data
             * stored across all ranks.
             */
            BufferedTable<contig_row_t> m_global;
            store_row_t m_source_row;

            typedef std::function<void(const store_row_t &, contig_row_t &)> loading_fn_t;
            const loading_fn_t m_loading_fn;

            PartSharedRowSet(const Communicator &comm, str_t name, contig_row_t contig_row, loading_fn_t loading_fn) :
                    DynamicRowSet(comm, name),
                    m_local("Dynamic shared row set \"" + name + "\" (local)", contig_row),
                    m_global("Dynamic shared set \"" + name + "\" (global)", contig_row),
                    m_source_row(comm.m_store.m_row), m_loading_fn(loading_fn) {
                m_local.push_back(1);
                m_global.push_back(1);
            }

            void update_data() {
                /*
                 * assumes that m_counts and m_displs have already been filled.
                 *
                 * The elements of m_irows determine the indices of the dynamic rows stored on this rank.
                 * First, we use these indices to randomly access the source MappedTable, and copy rows
                 * into the local contiguous table (m_local)
                 */
                DEBUG_ASSERT_TRUE_ALL(nrow(),
                                      logging::format(
                                              "Total number of rows across all ranks should be non-zero (\"{}\").",
                                              m_name));
                m_local.clear();
                auto nrow = m_displs.back() + m_counts.back();
                m_local.push_back(m_counts[mpi::irank()]);
                auto &local_row = m_local.m_row;
                local_row.restart();
                for (auto &irow: m_irows) {
                    m_source_row.jump(irow);
                    m_loading_fn(m_source_row, local_row);
                    local_row.step();
                }
                DEBUG_ASSERT_EQ(local_row.index(), m_counts[mpi::irank()], "not all local rows filled");

                m_global.clear();
                m_global.push_back(nrow);
                /*
                 * convert from units of rows to datawords...
                 */
                for (auto &i: m_counts) i *= m_local.row_size();
                for (auto &i: m_displs) i *= m_local.row_size();
                mpi::all_gatherv(m_local.begin(), m_counts[mpi::irank()], m_global.begin(), m_counts, m_displs);
                /*
                 * ... and back again
                 */
                for (auto &i: m_counts) i /= m_local.row_size();
                for (auto &i: m_displs) i /= m_local.row_size();
            }

            void update() override {
                DynamicRowSet::update();
                update_data();
            }
        };

        /**
         * A special case of PartSharedRowSet where all row data is copied, not just a subset
         */
        struct SharedRowSet : public PartSharedRowSet<store_row_t> {

            typedef typename PartSharedRowSet<store_row_t>::loading_fn_t loading_fn_t;

            static loading_fn_t make_copy_row_fn() {
                return [](const store_row_t &source, store_row_t &local) { local.copy_in(source); };
            }

            SharedRowSet(const Communicator &comm, str_t name) :
                    PartSharedRowSet<store_row_t>(comm, name, comm.m_store.m_row, make_copy_row_fn()) {}
        };

        struct SharedRow : public SharedRowSet {
            using DynamicRowSet::m_ra;
            using DynamicRowSet::m_source;
            using DynamicRowSet::m_irows;
            using DynamicRowSet::clear;
            using DynamicRowSet::add_;
            using DynamicRowSet::nrow_;
            using DynamicRowSet::nrow;
            using DynamicRowSet::update;
            using DynamicRowSet::m_name;
            using DynamicRowSet::m_ranks_with_any_rows;
            using SharedRowSet::m_global;

            uint_t m_iblock_ra;

            SharedRow(const Communicator &comm, TableBase::Loc loc, str_t name) :
                    SharedRowSet(comm, name) {
                redefine(loc);
                DEBUG_ASSERT_EQ(loc.is_mine(), comm.m_store.is_protected(),
                                "the SharedRow's table should be protected");
            }

            void redefine(TableBase::Loc newloc) {
                clear();
                m_global.m_row.restart();
                auto row = m_source.m_store.m_row;
                if (newloc.is_mine()) {
                    add_(newloc.m_irow);
                    row.jump(newloc.m_irow);
                    m_iblock_ra = m_source.m_ra.get_block(row);
                }
                mpi::bcast(m_iblock_ra, newloc.m_irank);
                update();
                logging::debug("Shared row \"{}\" is now in block {} on rank {}", m_name, m_iblock_ra, newloc.m_irank);
            }

            uint_t irank() const {
                return m_ranks_with_any_rows.size() ? m_ranks_with_any_rows[0] : ~0ul;
            }

            void update() override {
                auto irank_initial = irank();
                SharedRowSet::update();
                DEBUG_ASSERT_EQ(nrow(), 1ul, "Total number of rows across all ranks should be 1.");
                REQUIRE_EQ_ALL(m_ranks_with_any_rows.size(), 1ul, "Only one rank should have a row");
                auto irank_final = irank();
                if (irank_initial == ~0ul)
                    logging::info("Dynamic row \"{}\" is in block {} of {}, which is initially stored on rank {}",
                                  m_name, m_iblock_ra, this->m_ra.m_nblock, irank_final);
                else if (irank_initial != irank_final)
                    logging::info("Dynamic row \"{}\" moved from rank {} to rank {}",
                                  m_name, irank_initial, irank_final);
            }

            bool is_mine() const {
                ASSERT(m_ranks_with_any_rows.size() == 1);
                ASSERT(m_ranks_with_any_rows[0] == mpi::irank());
                return DynamicRowSet::nrow_();
            }

            bool is_same(const store_row_t &row) const {
                return is_mine() ? row.index() == *m_irows.cbegin() : false;
            }
        };

        /**
         * @param name
         *  base name for all buffers created here
         * @param store_nrow_est
         *  estimate of the number of rows that will ultimately be required in the store table
         * @param store_exp_fac
         *  fractional over-allocation to make when store buffer is resized
         * @param comm_nrow_est
         *  estimate of the number of rows per rank that will ultimately be required in the send table of the CommunicatingPair
         * @param comm_exp_fac
         *  fractional over-allocation to make when CommunicatingPair buffers are resized
         * @param store
         *  store table instance
         * @param send
         *  send table instance
         * @param nblock_ra
         *  number of RankAllocator blocks per rank
         * @param period_ra
         *  number of cycles between rank reallocation attempts
         * @param acceptable_imbalance
         *  fractional imbalance of work between busiest and laziest ranks acceptable to the rank allocator
         * @param nnull_updates_deactivate
         *  number of consecutive acceptably-imbalanced periods required for dynamic load balancing to be deactivated
         */
        Communicator(str_t name, uint_t store_nrow_est, double store_exp_fac, uint_t comm_nrow_est,
                     double comm_exp_fac, const store_table_t &store, const send_table_t &send, uint_t nblock_ra,
                     uint_t period_ra, double acceptable_imbalance, uint_t nnull_updates_deactivate) :
                m_store(name + " store", store),
                m_comm(name, comm_nrow_est, comm_exp_fac, send),
                m_lb(name, m_store, nblock_ra, period_ra, acceptable_imbalance, nnull_updates_deactivate),
                m_name(name) {
            logging::info("Initially allocating {} per rank for store buffer of \"{}\" ",
                          string::memsize(store_nrow_est * m_store.row_size()), name);
            m_store.resize(store_nrow_est, 0.0);
            m_store.set_expansion_factor(store_exp_fac);
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
#if 0
        Communicator(str_t name, uint_t store_nrow_crude_est, uint_t comm_nrow_crude_est,
                     const conf::Buffers &buf_opts, const conf::LoadBalancing &ra_opts,
                     const store_table_t &store, const send_table_t &send) :
                Communicator(name, std::max(1ul, store_nrow_crude_est) * buf_opts.m_store_fac_init, buf_opts.m_store_exp_fac,
                             std::max(1ul, comm_nrow_crude_est) * buf_opts.m_comm_fac_init, buf_opts.m_comm_exp_fac,
                             store, send,
                             ra_opts.m_nblock_per_rank * mpi::nrank(), ra_opts.m_period, ra_opts.m_acceptable_imbalance,
                             ra_opts.m_nnull_updates_deactivate) {}
#endif

        virtual ~Communicator() {}

        /**
         * @param key
         *  m_store-mapping key instance
         * @return
         *  rank index assigned to the key specified
         */
        uint_t get_rank(const key_field_t &key) const {
            return m_ra.get_rank(key);
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
#endif
}


#endif //M7_REDISTRIBUTOR_H