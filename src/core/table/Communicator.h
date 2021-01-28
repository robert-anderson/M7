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


template <typename table_t>
class CommunicatingPair {

    BufferedTableArray<table_t> m_send;
    BufferedTable<table_t> m_recv;
    double m_buffer_expansion_factor;

public:
    template<typename ...Args>
    CommunicatingPair(std::string name, double buffer_expansion_factor, const table_t& table):
            m_send(name+" send", mpi::nrank(), table), m_recv(name+" recv", table),
            m_buffer_expansion_factor(buffer_expansion_factor){
        m_send.set_expansion_factor(m_buffer_expansion_factor);
        m_recv.set_expansion_factor(m_buffer_expansion_factor);
    }

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
        m_recv.resize(nrow*mpi::nrank());
    }

    void expand(size_t nrow, double expansion_factor) {
        m_send.expand(nrow, expansion_factor);
        m_recv.expand(nrow*mpi::nrank(), expansion_factor);
    }

    void expand(size_t nrow) {
        m_send.expand(nrow);
        m_recv.expand(nrow*mpi::nrank());
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

        if (recv_dsize > static_cast<const Table&>(recv()).bw_dsize()){
            /*
             * the recv table is full
             * this expansion by a factor is done explicitly here, because we
             * want to expand relative to size of the incoming data, not the
             * current size of the buffer
             */
            m_recv.resize(std::ceil((1.0+m_buffer_expansion_factor)*recv_nrow));

        }

        if (!m_send.dbegin()) mpi::stop_all("Send buffer is not allocated!");
        if (!m_recv.dbegin()) mpi::stop_all("Recv buffer is not allocated!");

        auto tmp = mpi::all_to_allv(m_send.dbegin(), sendcounts, senddispls,
                                    m_recv.dbegin(), recvcounts, recvdispls);
        /*
         * check that the data addressed to this rank from this rank has been copied correctly
         */
        ASSERT(!send(mpi::irank()).dbegin() or std::memcmp(
                (void*) (send(mpi::irank()).dbegin()),
                (void*) (recv().dbegin()+recvdispls[mpi::irank()]),
                recvcounts[mpi::irank()]*defs::nbyte_data)==0);

        if (!tmp) throw std::runtime_error("MPI AllToAllV failed");

        recv().m_hwm = recv_nrow;
        m_send.clear();
    }

    void set_expansion_factor(double f){
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

    typedef typename store_t::base_table_t table_t;
    typedef typename store_t::key_field_t key_field_t;
    typedef typename key_field_t::view_t view_t;
    mutable RankAllocator<key_field_t> m_ra;
    BufferedTable<store_t> m_store;
    CommunicatingPair<comm_t> m_comm;
    double m_buffer_expansion_factor;

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
         * map row index in source -> row index in m_local.
         *
         * We could have achieved this mapping by extending the table_t type with an additional
         * field of type fields::Number<size_t>, then extending MappedTable with the resultant
         * type as a template argument. This has the downside that rows cannot be copied wholesale
         * from m_source into/out of such a table, since the rows would be longer in the extended
         * table_t. Better to use the stl-provided map to point into a table of the same format.
         */
        std::map<size_t, size_t> m_map;
        /*
         * once gathered, the major ordering of m_all will be MPI rank index ascending. m_counts stores
         * the number of rows tracked on each rank and m_displs stores the row index of m_all where
         * each rank's rows begin.
         */
        defs::inds m_counts;
        defs::inds m_displs;
        /*
         * keep track of which ranks have any rows
         */
        bool m_have_any_rows;
        defs::inds m_ranks_with_any_rows;

    public:
        DynamicRowSet(const Communicator &comm) : ra_t::Dynamic(comm.m_ra),
                m_source(comm.m_store),
                m_lc("Rank-local dynamic row set rows", static_cast<const table_t &>(comm.m_store)),
                m_ac("All dynamic row set rows", static_cast<const table_t &>(comm.m_store)),
                m_counts(mpi::nrank(), 0ul),
                m_displs(mpi::nrank(), 0ul) {
            m_ranks_with_any_rows.reserve(mpi::nrank());
        }

        size_t nrow() const {
            return m_map.size();
        }

        void add(size_t irow) {
            m_map[irow] = nrow();
        }

    private:
        void populate_local() {
            static_cast<Table &>(m_lc).clear();
            m_lc.push_back(nrow());
            /*
             * local row indices must be updated in case we sent a block
             */
            size_t irow_local = 0ul;
            for (auto &pair : m_map) {
                pair.second = irow_local++;
                static_cast<Table &>(m_lc).copy_row_in(m_source, pair.first, pair.second);
            }
        }

        void gatherv() {
            m_ac.clear();
            mpi::all_gather(nrow(), m_counts);
            mpi::counts_to_displs_consec(m_counts, m_displs);
            auto nrow = m_displs.back() + m_counts.back();
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
            m_have_any_rows = m_counts[mpi::irank()];
            m_ranks_with_any_rows.clear();
            for (auto &i : m_counts) {
                if (i) m_ranks_with_any_rows.push_back(mpi::irank());
            }
        }

    public:
        virtual void update() {
            populate_local();
            gatherv();
        }

        void on_row_send_(size_t irow) override {
            m_map.erase(irow);
        }

        void on_row_recv_(size_t irow) override {
            ASSERT(m_map.find(irow) == m_map.end());
            m_map[irow] = nrow();
        }

    };

    DynamicRowSet dynamic_row_set(RankAllocator<key_field_t> &ra) {
        return DynamicRowSet(ra, *this);
    }


    struct DynamicRow : public DynamicRowSet {
        using typename DynamicRowSet::ra_t;
        using ra_t::Dynamic::m_ra;
        using DynamicRowSet::m_source;
        using DynamicRowSet::m_map;
        using DynamicRowSet::add;
        using DynamicRowSet::nrow;
        using DynamicRowSet::update;
        using DynamicRowSet::m_have_any_rows;
        using DynamicRowSet::m_ranks_with_any_rows;

        std::string m_name;

        DynamicRow(const Communicator &comm, Table::Loc loc, std::string name) :
                DynamicRowSet(comm), m_name(std::move(name)) {
            if (loc.is_mine()) add(loc.m_irow);
            update();
            if (m_ranks_with_any_rows.size()!=1) mpi::stop_all("Only one rank should have a row");
        }

        Table::Loc location() const {
            if (m_ranks_with_any_rows.empty()) return {~0ul, ~0ul};
            return {m_ranks_with_any_rows[0], m_map.begin()->first};
        }

        void update() override {
            auto initial_loc = location();
            DynamicRowSet::update();
            auto final_loc = location();
            if (!initial_loc) log::info("Dynamic row \"{}\" is initially stored on rank {}",
                                        m_name, final_loc.m_irank);
            else if (initial_loc != final_loc) log::info("Dynamic row \"{}\" moved from rank {} to rank {}",
                                                    m_name, initial_loc.m_irank, final_loc.m_irank);
        }

        void change(Table::Loc loc) {
            if (loc.is_mine()) {
                ASSERT(m_ra.get_rank(m_source.m_key_field(loc.m_irow))==mpi::irank());
                m_map.clear();
                m_map = {loc.m_irow, 0ul};
            }
            update();
        }

        bool is_mine() const {
            ASSERT(m_have_any_rows);
            ASSERT(m_ranks_with_any_rows.size()==1);
            return DynamicRowSet::nrow();
        }
    };



    template<typename ...Args>
    Communicator(std::string name, double buffer_expansion_factor,
                 size_t nblock_ra, size_t period_ra,
                 const store_t &store, const comm_t &comm):
            m_ra(nblock_ra, period_ra),
            m_store(name + " store", store),
            m_comm(name, buffer_expansion_factor, comm),
            m_buffer_expansion_factor(buffer_expansion_factor) {
    }

    size_t get_rank(const view_t& key) const{
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
