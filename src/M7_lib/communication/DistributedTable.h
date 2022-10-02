//
// Created by rja on 01/10/22.
//

#ifndef M7_DISTRIBUTEDTABLE_H
#define M7_DISTRIBUTEDTABLE_H

#include "M7_lib/table/MappedTable.h"
#include "Distribution.h"
#include "SendRecv.h"

class DistribDependentBase {
    virtual void on_redistribution() = 0;
};

struct DistribBase {
    mutable std::set<DistribDependentBase*> m_dependents;
    ~DistribBase() {
        REQUIRE_TRUE_ALL(m_dependents.empty(), "Communicator should never be outlived by its dependents");
    }
};

struct DistribOptions {
    const uint_t m_nblock_per_rank;
    const uint_t m_period;
    const double m_acceptable_imbalance;
    const uint_t m_nnull_updates_deactivate;
};

template<typename row_t>
struct DistributedTable : MappedTable<row_t>, DistribBase {
    /**
     * the key field of the store table determines the block, and therefore MPI rank index, to which the row belongs
     */
    typedef typename KeyField<row_t>::type key_field_t;
    /**
     * current allocation of load balancing blocks to MPI rank indices
     */
    Distribution m_dist;
    /**
     * redistribution of rows is crucial to achieve load balancing. this pair must always have the same Row type as the
     * store table though
     */
    typedef send_recv::BasicSend<row_t> redist_pair_t;
    redist_pair_t m_redist;
    /**
     * row protection information must also be shared
     */
    typedef SingleFieldRow<field::Number<uint_t>> prot_level_row_t;
    typedef send_recv::BasicSend<prot_level_row_t> prot_level_pair_t;
//    prot_level_pair_t m_prot_level;
    /**
     * redistribution requires information about the amount of work done (by an arbitrary measure) by each block.
     */
    v_t<double> m_block_work_figures;

    static constexpr Sizing c_redist_sizing{1000ul, 1.0};

//    DistributedTable(str_t name, const row_t &row, uint_t nblock_per_rank, uint_t nbucket = 0ul,
//                     uint_t remap_nlookup = 0ul, double remap_ratio = 0.0) :
//            MappedTable<row_t>(row, nbucket, remap_nlookup, remap_ratio),
//            m_dist(nblock_per_rank*mpi::nrank()),
//            m_redist(name+" redistributor", row, c_redist_sizing),
//            m_prot_level(name+" protection level", {{}}, c_redist_sizing),
//            m_block_work_figures(m_dist.nblock()) {
//    }
    DistributedTable(const row_t &row, uint_t nblock_per_rank, uint_t nbucket = 0ul,
                     uint_t remap_nlookup = 0ul, double remap_ratio = 0.0) :
            MappedTable<row_t>(row, nbucket, remap_nlookup, remap_ratio),
            m_dist(nblock_per_rank*mpi::nrank()),
            m_redist("", row, c_redist_sizing),
//            m_prot_level("", SingleFieldRow<field::Number<uint_t>>(), c_redist_sizing),
            m_block_work_figures(m_dist.nblock()) {
    }

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

        /*
         * get a copy of the table<row_t>::m_row for traversal (cheap compared to loop)
         */
        auto row = this->row();
        {
            /*
             * loop over all rows,
             *  write those which no longer belong here into the communicating send table,
             *  clear them from the store
             */
            for (row.restart(); row.in_range(); row.step()) {
                if (row.key_field().is_zero()) continue;
                auto irank_owner = irank(row.key_field());
                if (!mpi::i_am(irank_owner)) {
                    m_redist.send(irank_owner).m_row.push_back_jump();
                    m_redist.send(irank_owner).m_row = row;
//                    m_prot_level.send(irank_owner).m_row.push_back_jump();
//                    m_prot_level.send(irank_owner).m_row.m_field = row.protection_level();
                    while (row.is_protected()) row.release();
                    clear_row(row);
                }
            }
            DEBUG_ASSERT_FALSE(m_redist.send(mpi::irank()).m_hwm, "rank should not be sending to itself");
        }
        m_redist.communicate();
//        m_prot_level.communicate();
        {
            auto recv_row = m_redist.recv().cursor();
//            auto prot_level_row = m_prot_level.recv().cursor();
//            prot_level_row.restart();
            for (recv_row.restart(); recv_row.in_range(); recv_row.step()) {
                const auto irank_owner = irank(recv_row.key_field());
                DEBUG_ONLY(irank_owner);
                DEBUG_ASSERT_TRUE(mpi::i_am(irank_owner), "recv_row sent to wrong rank!");
                insert(recv_row, row);
//                for (uint_t ilevel=0ul; ilevel < uint_t(prot_level_row.m_field); ++ilevel) row.protect();
//                prot_level_row.step();
            }
        }
        clear_work_figures();
    }
};


namespace buffered {
    template <typename row_t>
    struct DistributedTable : BufferedTable<row_t, ::DistributedTable<row_t>> {
        DistributedTable(const row_t &row, uint_t nblock_per_rank, uint_t nbucket = 0ul,
                         uint_t remap_nlookup = 0ul, double remap_ratio = 0.0):
            BufferedTable<row_t, ::DistributedTable<row_t>>(
                    ::DistributedTable<row_t>(row, nblock_per_rank, nbucket, remap_nlookup, remap_ratio)){}
    };
}


class DistribDependent : DistribDependentBase {
    const DistribBase& m_dist_base;
public:
    DistribDependent(const DistribBase& dist_base): m_dist_base(dist_base){
        m_dist_base.m_dependents.insert(this);
    }
    ~DistribDependent(){
        m_dist_base.m_dependents.erase(this);
    }
};



#endif //M7_DISTRIBUTEDTABLE_H
