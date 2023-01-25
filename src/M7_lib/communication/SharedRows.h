//
// Created by anderson on 27/07/2022.
//

#ifndef M7_SHAREDROWS_H
#define M7_SHAREDROWS_H

#include "Communicator.h"
#include "M7_lib/wavefunction/WalkerTable.h"

namespace shared_rows {
    /**
     * A set of rows which is copied to all ranks, responds correctly to rank reallocation and is protected from erasure
     */
    template<typename row_t>
    struct Set : public DistribDependent {
    protected:
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
         * contiguous send table for gathering into m_all
         */
        buffered::Table<row_t> m_gather_send;

        /**
         * row object for traversing the source storage table
         */
        row_t m_src_row;
        row_t m_send_row;

    public:
        Set(str_t name, const src_t& src, uintv_t irows = {}) :
                DistribDependent(src), m_src(src),
                m_name(name), m_all(name+" all rows", m_src.m_row, true),
                m_gather_send(name+" gather send", m_src.m_row),
                m_src_row(m_src.m_row), m_send_row(m_gather_send.m_row) {
            for (auto irow: irows) add_(irow);
            full_update();
            m_all.m_row.restart();
        }

        virtual ~Set() {}

    private:
        void prepare_gather() {
            if (m_gather_send.nrecord() < nrec_()) m_gather_send.resize(nrec_());
            m_gather_send.clear();
            for (auto irec: m_irecs) {
                m_src_row.jump(irec);
                m_send_row.push_back_jump();
                m_send_row.copy_in(m_src_row);
            }
        }

    public:

        /**
         * bring all rows into m_all
         * direct buffer copy, no rows have changed locations so no need to treat m_all as a mapped table
         */
        void update() {
            prepare_gather();
            m_all.TableBase::all_gatherv(m_gather_send);
        }

        void row_index_update() {
            REQUIRE_EQ_ALL(m_all.nrow_in_use(), m_irecs.size(),
                           "this kind of refresh is only possible when all rows have already been gathered");
            m_irecs = {};
            auto& row = m_all.m_row;
            for (row.restart(); row; ++row) {
                auto result = m_src.lookup(row.key_field(), m_src_row);
                if (result) {
                    const auto irec = static_cast<const Row&>(m_src_row).index();
                    m_irecs.insert(irec);
                }
            }
        }

        /**
         * include possible changes in basis order in the update
         */
        void full_update() {
            prepare_gather();
            m_all.all_gatherv(m_gather_send);
        }

        uint_t nrec_() const {
            return m_irecs.size();
        }

//        uint_t nrec() const {
//            return mpi::all_sum(nrow_());
//        }

        void clear() {
            for (const auto& irow: m_irecs) m_src.unprotect(irow);
            m_irecs.clear();
        }

        void add_(uint_t irow) {
            m_src.protect(irow);
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
        using Set<row_t>::m_irecs;

        Single(str_t name, const src_t& src, TableBase::Loc loc): Set<row_t>(name, src, loc) {}

        void redefine(TableBase::Loc loc) {
            clear();
            if (loc.is_mine()) add_(loc.m_irec);
            full_update();
        }

        uint_t irec() const {
            return m_irecs.empty() ? ~0ul : *m_irecs.cbegin();
        }

        bool is_mine() const {
            return irec() != ~0ul;
        }
    };

    struct Walker : Single<::Walker> {

        Walker(str_t name, const src_t& src, TableBase::Loc loc): Single<::Walker>(name, src, loc) {}

        const field::Mbf &mbf() const {
            return m_all.m_row.m_mbf;
        }

        const wf_t &weight(uint_t ipart) const {
            return m_all.m_row.m_weight[ipart];
        }

        /**
         * this method includes the current weight in the average, bringing the normalized average up to date.
         * @param icycle
         *  cycle on which average is being used
         * @param ipart
         *  wf part index
         * @return
         *  normalized average weight
         */
        wf_t norm_average_weight(uint_t icycle, uint_t ipart) const {
            auto unnorm = m_all.m_row.m_average_weight[ipart]+m_all.m_row.m_weight[ipart];
            return unnorm/static_cast<wf_comp_t>(m_all.m_row.occupied_ncycle(icycle));
        }
    };
}

#endif //M7_SHAREDROWS_H
