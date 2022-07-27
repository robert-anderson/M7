//
// Created by anderson on 27/07/2022.
//

#ifndef M7_SHAREDROWS_H
#define M7_SHAREDROWS_H

#include "Communicator.h"

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

#endif //M7_SHAREDROWS_H
