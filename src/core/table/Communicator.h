//
// Created by rja on 09/11/2020.
//

#ifndef M7_COMMUNICATOR_H
#define M7_COMMUNICATOR_H

#include "BufferedTable.h"
#include "BufferedTableArray.h"
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/io/Logging.h"
#include <set>
#include <src/core/parallel/RankAllocator.h>
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
template<typename row_t, bool mapped=false>
class CommunicatingPair {
public:
    typedef BufferedTableArray<row_t, mapped> send_t;
    typedef BufferedTable<row_t> recv_t;
    typedef typename send_t::table_t send_table_t;
    typedef typename recv_t::table_t recv_table_t;
private:
    send_t m_send;
    recv_t m_recv;
    double m_buffer_expansion_factor;

public:
    /*
     * number of rows sent and recvd in the last call to communicate()
     * retained for stats reasons
     */
    defs::inds m_last_send_counts;
    size_t m_last_recv_count = 0ul;
    template<typename ...Args>
    CommunicatingPair(std::string name, double buffer_expansion_factor, const send_table_t &send):
            m_send(name + " send", mpi::nrank(), send),
            m_recv(name + " recv", recv_table_t(send.m_row)),
            m_buffer_expansion_factor(buffer_expansion_factor) {
        m_send.set_expansion_factor(m_buffer_expansion_factor);
        m_recv.set_expansion_factor(m_buffer_expansion_factor);
    }

    size_t row_dsize() const {
        return static_cast<const TableBase &>(m_recv).m_row_dsize;
    }

    send_t &send() {
        return m_send;
    }

    const send_t &send() const {
        return m_send;
    }

    typename send_t::table_t &send(const size_t &i) {
        return m_send[i];
    }

    const typename send_t::table_t &send(const size_t &i) const {
        return m_send[i];
    }

    Table<row_t> &recv() {
        return m_recv;
    }

    const Table<row_t> &recv() const {
        return m_recv;
    }

    void resize(size_t nrow) {
        m_send.resize(nrow);
        m_recv.resize(nrow * mpi::nrank());
    }

    void expand(size_t nrow, double expansion_factor) {
        m_send.expand(nrow, expansion_factor);
        m_recv.expand(nrow * mpi::nrank(), expansion_factor);
    }

    void expand(size_t nrow) {
        m_send.expand(nrow);
        m_recv.expand(nrow * mpi::nrank());
    }


    void communicate() {
        m_last_send_counts = m_send.hwms();
        defs::inds sendcounts(m_last_send_counts);
        for (auto &it: sendcounts) it *= row_dsize();
        defs::inds recvcounts(mpi::nrank(), 0ul);

        log::debug_("Sending {} rows", utils::to_string(m_last_send_counts));

        mpi::all_to_all(sendcounts, recvcounts);

        auto senddispls = m_send.displs();
        defs::inds recvdispls(mpi::nrank(), 0ul);
        for (size_t i = 1ul; i < mpi::nrank(); ++i)
            recvdispls[i] = recvdispls[i - 1] + recvcounts[i - 1];
        auto recv_dsize = recvdispls.back() + recvcounts.back();
        m_last_recv_count = recv_dsize / row_dsize();

        if (recv_dsize > static_cast<const TableBase &>(recv()).bw_dsize()) {
            /*
             * the recv table is full
             * this expansion by a factor is done explicitly here, because we
             * want to expand relative to size of the incoming data, not the
             * current size of the buffer
             */
            m_recv.resize(std::ceil((1.0 + m_buffer_expansion_factor) * m_last_recv_count));
        }

        REQUIRE_TRUE_ALL(m_send.dbegin(), "Send buffer is not allocated on all ranks!");
        REQUIRE_TRUE_ALL(m_recv.dbegin(), "Recv buffer is not allocated on all ranks!");

        auto tmp = mpi::all_to_allv(m_send.dbegin(), sendcounts, senddispls,
                                    m_recv.dbegin(), recvcounts, recvdispls);
        /*
         * check that the data addressed to this rank from this rank has been copied correctly
         */
        ASSERT(!send(mpi::irank()).dbegin() or std::memcmp(
                (void *) (send(mpi::irank()).dbegin()),
                (void *) (recv().dbegin() + recvdispls[mpi::irank()]),
                recvcounts[mpi::irank()] * defs::nbyte_data) == 0);

        if (!tmp) throw std::runtime_error("MPI AllToAllV failed");

        recv().m_hwm = m_last_recv_count;
        m_send.clear();
    }

    void set_expansion_factor(double f) {
        m_send.set_expansion_factor(f);
        m_recv.set_expansion_factor(f);
    }
};

/**
 * Combines the CommunicatingPair with a persistent storage table and a RankAllocator
 * @tparam store_row_t
 *  The Row class-derived data layout of the storage table
 * @tparam comm_row_t
 *  The Row class-derived data layout of the tables within the CommunicatingPair
 * @tparam mapped_comm
 *  optional mapping of CommunicatingPair
 */
template<typename store_row_t, typename comm_row_t, bool mapped_comm=false>
struct Communicator {
    static_assert(std::is_base_of<Row, store_row_t>::value, "Template arg must be derived from Row");
    static_assert(std::is_base_of<Row, comm_row_t>::value, "Template arg must be derived from Row");

    typedef typename KeyField<store_row_t>::type key_field_t;
    typedef BufferedTable<store_row_t, true> store_t;
    typedef CommunicatingPair<comm_row_t, mapped_comm> comm_t;
    typedef typename store_t::table_t store_table_t;
    typedef typename comm_t::send_t::table_t send_table_t;

    store_t m_store;
    comm_t m_comm;
    mutable RankAllocator<store_row_t> m_ra;
    std::string m_name;
    double m_buffer_expansion_factor;

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
        std::set<size_t> m_irows;
        /**
         * the counts and displs determine a global index for a local dynamic row. The ith element of m_irows on MPI
         * rank j has global index m_displs[j]+i
         */
        defs::inds m_counts;
        defs::inds m_displs;
        defs::inds m_ranks_with_any_rows;

        std::string m_name;
        /*
         * row = row in mapped table
         * trow = row among transferred rows
         * drow = dynamic row index among transferred rows
         */
        const int m_ntrow_to_track_p2p_tag;
        const int m_itrows_to_track_p2p_tag;
        defs::inds m_idrows;
        size_t m_itrow;
        size_t m_ndrow_found;

    public:
        DynamicRowSet(const Communicator &comm, std::string name) :
                RankDynamic(comm.m_ra),
                RowProtector(comm.m_store),
                m_source(comm),
                m_counts(mpi::nrank(), 0ul),
                m_displs(mpi::nrank(), 0ul),
                m_name(name),
                m_ntrow_to_track_p2p_tag(mpi::new_p2p_tag()),
                m_itrows_to_track_p2p_tag(mpi::new_p2p_tag()) {
            log::debug("P2P tag for number of dynamic row indices to transfer for \"{}\": {}", m_name,
                       m_ntrow_to_track_p2p_tag);
            log::debug("P2P tag for array of dynamic row indices to transfer for \"{}\": {}", m_name,
                       m_itrows_to_track_p2p_tag);
            m_ranks_with_any_rows.reserve(mpi::nrank());
        }

        virtual ~DynamicRowSet(){}

        size_t nrow_() const {
            return m_irows.size();
        }

        size_t nrow() const {
            return mpi::all_sum(nrow_());
        }

        void clear() {
            for (const auto& irow: m_irows) RowProtector::release(irow);
            m_irows.clear();
        }

        void add_(size_t irow) {
            RowProtector::protect(irow);
            m_irows.insert(irow);
        }

        void add_(const store_row_t& row) {
            add_(row.index());
        }

    public:
        virtual void update() {
            mpi::all_gather(nrow_(), m_counts);
            mpi::counts_to_displs_consec(m_counts, m_displs);
            m_ranks_with_any_rows.clear();
            for (size_t irank=0; irank<mpi::nrank(); ++irank){
                if (m_counts[irank]) m_ranks_with_any_rows.push_back(irank);
            }
        }

        bool has_row(size_t irow) override {
            return m_irows.find(irow)!=m_irows.end();
        }

        void before_block_transfer(const defs::inds &irows_send, size_t irank_send, size_t irank_recv) override {
            size_t nrow_transfer;
            if (mpi::i_am(irank_send)) {
                m_idrows.clear();
                // look for Dynamic rows among those being transferred
                size_t itrow = 0ul;
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
                    log::debug_("No dynamic rows to transfer for \"{}\"", m_name);
                    return;
                }
                log::debug_("Sending {} dynamic rows for \"{}\" from rank {} to {}",
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
                    log::debug_("No dynamic rows to transfer for \"{}\"", m_name);
                    return;
                }
                log::debug_("Recving {} dynamic rows for \"{}\" to rank {} from {}",
                            nrow_transfer, m_name, irank_recv, irank_send);
                m_idrows.resize(nrow_transfer, ~0ul);
                mpi::recv(m_idrows.data(), nrow_transfer, irank_send, m_itrows_to_track_p2p_tag);
                log::debug_("idrows: {}", utils::to_string(m_idrows));
            }
        }

        void on_row_recv_(size_t irow) override {
            // check if there's any more rows we may need to track
            if (m_ndrow_found < m_idrows.size()) {
                ASSERT(m_ndrow_found < m_idrows.size());
                const auto next_itrow = m_idrows[m_ndrow_found];
                if (m_itrow == next_itrow) {
                    // this transferred row is dynamic
                    DEBUG_ASSERT_TRUE(m_irows.find(irow)==m_irows.end(), "Transferred dynamic row shouldn't already be here!");
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

        typedef std::function<void(const store_row_t&, contig_row_t&)> loading_fn_t;
        const loading_fn_t m_loading_fn;

        PartSharedRowSet(const Communicator &comm, std::string name, contig_row_t contig_row, loading_fn_t loading_fn) :
                DynamicRowSet(comm, name),
                m_local("Dynamic shared row set \"" + name + "\" (local)", contig_row),
                m_global("Dynamic shared set \"" + name + "\" (global)", contig_row),
                m_source_row(comm.m_store.m_row), m_loading_fn(loading_fn){
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
            DEBUG_ASSERT_TRUE_ALL(nrow(), log::format("Total number of rows across all ranks should be non-zero (\"{}\").", m_name));
            m_local.clear();
            auto nrow = m_displs.back() + m_counts.back();
            m_local.push_back(m_counts[mpi::irank()]);
            auto& local_row = m_local.m_row;
            local_row.restart();
            for (auto &irow : m_irows){
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
            for (auto &i : m_counts) i *= m_local.m_row_dsize;
            for (auto &i : m_displs) i *= m_local.m_row_dsize;
            mpi::all_gatherv(m_local.dbegin(), m_counts[mpi::irank()], m_global.dbegin(), m_counts, m_displs);
            /*
             * ... and back again
             */
            for (auto &i : m_counts) i /= m_local.m_row_dsize;
            for (auto &i : m_displs) i /= m_local.m_row_dsize;
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
            return [](const store_row_t& source, store_row_t& local){local.copy_in(source);};
        }

        SharedRowSet(const Communicator &comm, std::string name):
                PartSharedRowSet<store_row_t>(comm, name, comm.m_store.m_row, make_copy_row_fn()){}
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

        size_t m_iblock_ra;

        SharedRow(const Communicator &comm, TableBase::Loc loc, std::string name) :
                SharedRowSet(comm, name) {
            redefine(loc);
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
            log::debug("Shared row \"{}\" is now in block {} on rank {}", m_name, m_iblock_ra, newloc.m_irank);
        }

        size_t irank() const {
            return m_ranks_with_any_rows.size()?m_ranks_with_any_rows[0]:~0ul;
        }

        void update() override {
            auto irank_initial = irank();
            SharedRowSet::update();
            DEBUG_ASSERT_EQ(nrow(), 1ul, "Total number of rows across all ranks should be 1.");
            REQUIRE_EQ_ALL(m_ranks_with_any_rows.size(), 1ul, "Only one rank should have a row");
            auto irank_final = irank();
            if (irank_initial == ~0ul)
                log::info("Dynamic row \"{}\" is in block {} of {}, which is initially stored on rank {}",
                          m_name, m_iblock_ra, this->m_ra.m_nblock, irank_final);
            else if (irank_initial != irank_final)
                log::info("Dynamic row \"{}\" moved from rank {} to rank {}",
                          m_name, irank_initial, irank_final);
        }

        bool is_mine() const {
            ASSERT(m_ranks_with_any_rows.size() == 1);
            ASSERT(m_ranks_with_any_rows[0]==mpi::irank());
            return DynamicRowSet::nrow_();
        }

        bool is_same(const store_row_t& row) const {
            return is_mine() ? row.index() == *m_irows.cbegin() : false;
        }
    };


    Communicator(std::string name, double buffer_expansion_factor,
                 size_t nblock_ra, size_t period_ra,
                 const store_table_t &store, const send_table_t &send,
                 double acceptable_imbalance):
            m_store(name + " store", store),
            m_comm(name, buffer_expansion_factor, send),
            m_ra(m_store, nblock_ra, period_ra, acceptable_imbalance),
            m_name(name),
            m_buffer_expansion_factor(buffer_expansion_factor) {
    }

    virtual ~Communicator(){}

    size_t get_rank(const key_field_t &key) const {
        return m_ra.get_rank(key);
    }

    typename comm_t::send_t &send() {
        return m_comm.send();
    }

    const typename comm_t::send_t &send() const {
        return m_comm.send();
    }

    send_table_t &send(const size_t &i) {
        return m_comm.send(i);
    }

    const send_table_t &send(const size_t &i) const {
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

    void set_expansion_factor(double f) {
        m_buffer_expansion_factor = f;
        m_comm.set_expansion_factor(f);
    }

};

#endif //M7_COMMUNICATOR_H
