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
    template<typename row_t>
    struct Set : public DistribDependent {
        /**
         * the table with the definitive record values to be copied and sent to all ranks
         */
        typedef buffered::DistributedTable<row_t> src_t;
        const src_t& m_src;
        /**
         * set of record indices stored on this rank
         */
        std::set<uint_t> m_irecs;
        /**
         * name to use for this row set in detailed logging
         */
        str_t m_name;
        /**
         * all rows gathered
         */
        buffered::MappedTable<row_t> m_all;
        /**
         * send / recv tables for gathering into m_all
         */
        buffered::Table<row_t> m_gather_send;
        buffered::Table<row_t> m_gather_recv;

        /**
         * row object for traversing the source storage table
         */
        row_t m_src_row;
        row_t m_send_row;

    public:
        Set(str_t name, const src_t& src, uintv_t irecs = {}) :
                DistribDependent(src), m_src(src),
                m_name(name), m_all(name+" all rows", m_src.m_row),
                m_gather_send(name+" gather send", m_src.m_row),
                m_gather_recv(name+" gather recv", m_src.m_row),
                m_src_row(m_src.m_row), m_send_row(m_gather_send.m_row) {
            m_all.m_row.restart();
            for (auto irow: irecs) add_(irow);
            full_update();
        }

        virtual ~Set() {}

        /**
         * bring all rows into m_gather_recv
         */
        void all_gatherv() {
            m_gather_send.clear();
            for (auto irec: m_irecs) {
                m_src_row.jump(irec);
                m_send_row.push_back_jump();
                m_send_row.copy_in(m_src_row);
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
            REQUIRE_EQ_ALL(m_all.m_hwm, m_irecs.size(),
                           "this kind of refresh is only possible when all rows have already been gathered");
            m_irecs = {};
            auto& row = m_all.m_row;
            for (row.restart(); row.in_range(); row.step()) {
                auto result = m_src.lookup(row.key_field(), m_src_row);
                if (result) {
                    const auto irec = static_cast<const Row&>(m_src_row).index();
                    m_irecs.insert(irec);
                }
            }
        }

        void full_update() {
            all_gatherv();
            m_all.clear();
            auto& recv_row = m_gather_recv.m_row;
            for (recv_row.restart(); recv_row.in_range(); recv_row.step()) m_all.insert(recv_row);
        }

        uint_t nrec_() const {
            return m_irecs.size();
        }

//        uint_t nrec() const {
//            return mpi::all_sum(nrow_());
//        }

        void clear() {
            for (const auto& irow: m_irecs) static_cast<const TableBase&>(m_src).release(irow);
            m_irecs.clear();
        }

        void add_(uint_t irow) {
            static_cast<const TableBase&>(m_src).protect(irow);
            m_irecs.insert(irow);
        }

        void add_(const row_t& row) {
            add_(row.index());
        }

    private:
        void on_redistribution() override {
            row_index_update();
            full_update();
        }
    };

    template<typename row_t>
    struct Single : Set<row_t> {
        typedef buffered::DistributedTable<row_t> src_t;
        using Set<row_t>::m_all;
        using Set<row_t>::clear;
        using Set<row_t>::add_;
        using Set<row_t>::full_update;

        Single(str_t name, const src_t& src, TableBase::Loc loc): Set<row_t>(name, src, loc) {}

        void redefine(TableBase::Loc loc) {
            clear();
            if (loc.is_mine()) add_(loc.m_irec);
            full_update();
        }
    };
}

#endif //M7_SHAREDROWS_H
