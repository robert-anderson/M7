//
// Created by rja on 09/11/2020.
//

#ifndef M7_COMMUNICATOR_H
#define M7_COMMUNICATOR_H

#include "BufferedTableArray.h"
#include "BufferedTable.h"
#include "src/core/parallel/MPIWrapper.h"
#include "src/core/io/Logging.h"
#include "MappedTable.h"
#include <set>


template<typename table_t>
class CommunicatingPair {

    BufferedTableArray<table_t> m_send;
    BufferedTable<table_t> m_recv;
    double m_buffer_expansion_factor;

public:
    template<typename ...Args>
    CommunicatingPair(std::string name, double buffer_expansion_factor, const table_t &table):
            m_send(name + " send", mpi::nrank(), table), m_recv(name + " recv", table),
            m_buffer_expansion_factor(buffer_expansion_factor) {
        m_send.set_expansion_factor(m_buffer_expansion_factor);
        m_recv.set_expansion_factor(m_buffer_expansion_factor);
    }

    size_t row_dsize() const {
        return static_cast<const Table &>(m_recv).m_row_dsize;
    }

    BufferedTableArray<table_t> &send() {
        return m_send;
    }

    const BufferedTableArray<table_t> &send() const {
        return m_send;
    }

    table_t &send(const size_t &i) {
        return static_cast<table_t &>(m_send[i]);
    }

    const table_t &send(const size_t &i) const {
        return static_cast<const table_t &>(m_send[i]);
    }

    table_t &recv() {
        return static_cast<table_t &>(m_recv);
    }

    const table_t &recv() const {
        return static_cast<const table_t &>(m_recv);
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
        auto hwms = m_send.hwms();
        defs::inds sendcounts(hwms);
        for (auto &it: sendcounts) it *= row_dsize();
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

        if (recv_dsize > static_cast<const Table &>(recv()).bw_dsize()) {
            /*
             * the recv table is full
             * this expansion by a factor is done explicitly here, because we
             * want to expand relative to size of the incoming data, not the
             * current size of the buffer
             */
            m_recv.resize(std::ceil((1.0 + m_buffer_expansion_factor) * recv_nrow));
        }

        MPI_REQUIRE_ALL(m_send.dbegin(), "Send buffer is not allocated!");
        MPI_REQUIRE_ALL(m_recv.dbegin(), "Recv buffer is not allocated!");

        log::info_(utils::to_string(sendcounts));
        log::info_(utils::to_string(senddispls));
        log::info_(utils::to_string(recvcounts));
        log::info_(utils::to_string(recvdispls));
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

        recv().m_hwm = recv_nrow;
        m_send.clear();
    }

    void set_expansion_factor(double f) {
        m_send.set_expansion_factor(f);
        m_recv.set_expansion_factor(f);
    }
};

template<typename store_t, typename comm_t>
struct Communicator {
    static_assert(std::is_base_of<MappedTableBase, store_t>::value,
                  "Template arg must be derived from MappedTableBase");
    static_assert(std::is_base_of<Table, comm_t>::value,
                  "Template arg must be derived from Table");

    typedef typename store_t::table_t table_t;
    typedef typename store_t::key_field_t key_field_t;
    typedef typename key_field_t::view_t view_t;
    BufferedTable<store_t> m_store;
    CommunicatingPair<comm_t> m_comm;
    mutable RankAllocator<key_field_t> m_ra;
    std::string m_name;
    double m_buffer_expansion_factor;

    /**
     * The behaviour of this class is intricate.
     *
     * The base RankAllocator<key_field_t>::Dynamic has three pure virtual methods:
     *
     * 1. void on_row_send_(size_t irow)
     * 2. void before_block_transfer(size_t iblock, size_t irank_send, size_t irank_recv)
     * 3. void on_row_recv_(size_t irow)
     *
     * As indicated by the trailing _, functions 1 and 3 are called only on the affected
     * ranks, function 2 on the other hand is called on all ranks.
     *
     * std::function instances which wrap these are termed "callbacks". This is done so
     * that lists containing calls to these three methods in the
     * RankAllocator<key_field_t>::Dynamic can be called in succession within any given
     * context for all elements of RankAllocator<key_field_t>::m_dependents.
     *
     * m_dependents is a list of pointers to RankAllocator<key_field_t>::Dynamic, since
     * there may in general be many objects which must be dynamic with respect to
     * key block transfers.
     * E.g. Reference determinant (single row derived class "DynamicRow"),
     *      Trial wavefunction connected space
     *
     * //////////////////////////////////////////////////////////////////////////////////
     *
     * As the source MappedTable is scanned for keys belonging to the outgoing block index
     * in RankAllocator<key_field_t>::update, each time an affected row is found, its index
     * is appended to a vector.
     *
     * The Table class contains methods to send and recv rows in a way that is completely
     * unaware of keyfields and rank allocation. All that is required is a list of rows to
     * send and the indices of the sending and recving ranks. Appropriate overrides are in
     * place to accommodate the additional structure of the MappedTable. Both of these P2P
     * communication handlers accept an optional argument, a list of callbacks which
     * determine how other objects react to this transfer
     *
     * In Table::send_rows_, each callback is invoked just before each row is freed.
     *
     * each instance of
     * RankAllocator<key_field_t>::Dynamic attached to the RankAllocator<key_field_t> has
     * its on_row_send_ virtual method called, in the case
     * of the DynamicRowSet, this method looks
     *
     *
     *
     */
    struct DynamicRowSet : RankAllocator<key_field_t>::Dynamic {
        typedef RankAllocator<key_field_t> ra_t;
        /*
         * the mapped table which stores the definitive row values
         */
        const store_t &m_source;
        /*
         * the unmapped table which loads copies of rows between m_source (arbitrary order, non-
         * contiguous) and m_all (contiguous). lc = "local, contiguous"
         */
        BufferedTable<table_t> m_lc;
        /*
         * the unmapped table which holds copies of all mapped rows. these copies are refreshed
         * with a call to refresh method. This table is intended for reading rows from all MPI ranks,
         * as such there is no machanism for committing changes to m_source from m_all. It is
         * to be treated as a read-only copy. ac = "all, contiguous"
         */
        BufferedTable<table_t> m_ac;
        /*
         * set of dynamic row indices stored on this rank
         */
        std::set<size_t> m_irows;
        /*
         * once gathered, the major ordering of m_all will be MPI rank index ascending. m_counts stores
         * the number of rows tracked on each rank and m_displs stores the row index of m_all where
         * each rank's rows begin.
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
                ra_t::Dynamic(comm.m_ra),
                m_source(comm.m_store),
                m_lc("Dynamic row set \"" + name + "\" (local)",
                     static_cast<const table_t &>(comm.m_store)),
                m_ac("Dynamic row set \"" + name + "\" (all)",
                     static_cast<const table_t &>(comm.m_store)),
                m_counts(mpi::nrank(), 0ul),
                m_displs(mpi::nrank(), 0ul),
                m_ranks_with_any_rows(mpi::nrank(), 0ul),
                m_name(name),
                m_ntrow_to_track_p2p_tag(mpi::new_p2p_tag()),
                m_itrows_to_track_p2p_tag(mpi::new_p2p_tag()) {
            log::debug("P2P tag for number of dynamic row indices to transfer for \"{}\": {}", m_name,
                       m_ntrow_to_track_p2p_tag);
            log::debug("P2P tag for array of dynamic row indices to transfer for \"{}\": {}", m_name,
                       m_itrows_to_track_p2p_tag);
            m_lc.resize(1);
            m_ac.resize(1);
        }

        size_t nrow_() const {
            return m_irows.size();
        }

        size_t nrow() const {
            return mpi::all_sum(nrow_());
        }

        void add_(size_t irow) {
            m_irows.insert(irow);
        }


    protected:
        void update_counts() {
            mpi::all_gather(nrow_(), m_counts);
            mpi::counts_to_displs_consec(m_counts, m_displs);
            m_ranks_with_any_rows.clear();
            for (size_t irank=0; irank<mpi::nrank(); ++irank){
                if (m_counts[irank]) m_ranks_with_any_rows.push_back(irank);
            }
        }

        void update_data() {
            /*
             * It is expected that m_counts and m_displs have already been filled.
             *
             * The elements of m_irows determine the indices of the dynamic rows stored on this rank.
             * First, we use these indices to randomly access the source MappedTable, and copy rows
             * into the local contiguous table (m_lc)
             */
            MPI_ASSERT_ALL(nrow(), "Total number of rows across all ranks should be non-zero.");
            size_t irow_local = 0ul;
            m_lc.push_back(nrow_());
            for (auto &irow : m_irows) static_cast<Table &>(m_lc).copy_row_in(m_source, irow, irow_local++);
            ASSERT(irow_local==nrow_());

            m_ac.clear();
            mpi::all_gather(nrow_(), m_counts);
            mpi::counts_to_displs_consec(m_counts, m_displs);
            auto nrow = m_displs.back() + m_counts.back();
            ASSERT(nrow==this->nrow());
            m_ac.push_back(nrow);
            /*
             * convert from units of rows to datawords...
             */
            for (auto &i : m_counts) i *= m_source.m_row_dsize;
            for (auto &i : m_displs) i *= m_source.m_row_dsize;
            mpi::all_gatherv(m_lc.dbegin(), m_counts[mpi::irank()], m_ac.dbegin(), m_counts, m_displs);
            /*
             * ... and back again
             */
            for (auto &i : m_counts) i /= m_source.m_row_dsize;
            for (auto &i : m_displs) i /= m_source.m_row_dsize;
        }


    public:
        virtual void update() {
            update_counts();
            update_data();
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
                    MPI_ASSERT(m_irows.find(irow) == m_irows.end(), "Transferred dynamic row shouldn't already be here!");
                    m_irows.insert(irow);
                    ++m_ndrow_found;
                }
            }
            ++m_itrow;
        }

        void after_block_transfer() override {
            update();
        }

    };

    struct DynamicRow : public DynamicRowSet {
        using typename DynamicRowSet::ra_t;
        using ra_t::Dynamic::m_ra;
        using DynamicRowSet::m_source;
        using DynamicRowSet::m_irows;
        using DynamicRowSet::add_;
        using DynamicRowSet::nrow_;
        using DynamicRowSet::nrow;
        using DynamicRowSet::update;
        using DynamicRowSet::m_name;
        using DynamicRowSet::m_ranks_with_any_rows;

        size_t m_iblock;

        DynamicRow(const Communicator &comm, Table::Loc loc, std::string name) :
                DynamicRowSet(comm, name) {
            if (loc.is_mine()) {
                add_(loc.m_irow);
                m_iblock = m_ra.get_block(m_source.m_key_field(loc.m_irow));
            }
            mpi::bcast(m_iblock, loc.m_irank);
            update();
        }

        size_t irank() const {
            return m_ranks_with_any_rows.size()?m_ranks_with_any_rows[0]:~0ul;
        }

        void update() override {
            DynamicRowSet::update_counts();
            auto irank_initial = irank();
            DynamicRowSet::update_data();
            MPI_ASSERT(nrow() == 1, "Total number of rows across all ranks should be 1.");
            MPI_REQUIRE_ALL(m_ranks_with_any_rows.size() == 1, "Only one rank should have a row");
            auto irank_final = irank();
            if (irank_initial == ~0ul)
                log::info("Dynamic row \"{}\" is in block {} of {}, which is initially stored on rank {}",
                          m_name, m_iblock, m_ra.m_nblock, irank_final);
            else if (irank_initial != irank_final)
                log::info("Dynamic row \"{}\" moved from rank {} to rank {}",
                          m_name, irank_initial, irank_final);
        }

        void change(Table::Loc loc) {
            if (loc.is_mine()) {
                ASSERT(m_ra.get_rank(m_source.m_key_field(loc.m_irow)) == mpi::irank());
                m_irows.clear();
                m_irows = {loc.m_irow};
            }
            update();
        }

        bool is_mine() const {
            ASSERT(m_ranks_with_any_rows.size() == 1);
            ASSERT(m_ranks_with_any_rows[0]==mpi::irank());
            return DynamicRowSet::nrow_();
        }
    };


    template<typename ...Args>
    Communicator(std::string name, double buffer_expansion_factor,
                 size_t nblock_ra, size_t period_ra,
                 const store_t &store, const comm_t &comm):
            m_store(name + " store", store),
            m_comm(name, buffer_expansion_factor, comm),
            m_ra(m_store, m_store.m_key_field, nblock_ra, period_ra),
            m_name(name),
            m_buffer_expansion_factor(buffer_expansion_factor) {
    }

    size_t get_rank(const view_t &key) const {
        return m_ra.get_rank(key);
    }

    BufferedTableArray<comm_t> &send() {
        return m_comm.send();
    }

    const BufferedTableArray<comm_t> &send() const {
        return m_comm.send();
    }

    comm_t &send(const size_t &i) {
        return static_cast<comm_t &>(m_comm.send(i));
    }

    const comm_t &send(const size_t &i) const {
        return static_cast<const comm_t &>(m_comm.send(i));
    }

    comm_t &recv() {
        return static_cast<comm_t &>(m_comm.recv());
    }

    const comm_t &recv() const {
        return static_cast<const comm_t &>(m_comm.recv());
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
