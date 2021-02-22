//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include "src/core/field/FieldBase.h"
#include "src/core/table/MappedTable.h"
#include "src/core/parallel/Gatherable.h"
#include "src/core/io/Logging.h"
#include "MPIWrapper.h"
#include "Epoch.h"
#include <forward_list>
#include <algorithm>
#include <src/core/util/Timer.h>

struct RankAllocatorBase {
    struct Dependent {
        typedef RankAllocatorBase ra_t;
        ra_t& m_ra;
        typename std::list<Dependent*>::iterator m_it;
        Dependent(ra_t& ra): m_ra(ra), m_it(m_ra.add_dependent(this)) {}
        ~Dependent() {
            m_ra.erase_dependent(this);
        }

        virtual bool has_row(size_t irow) = 0;
        // called before block transfer is performed
        virtual void before_block_transfer(const defs::inds& irows_send, size_t irank_send, size_t irank_recv) = 0;
        // called after row is inserted into the MappedTable on the recving rank
        virtual void on_row_recv_(size_t irow) = 0;
        // called after block transfer is performed
        virtual void after_block_transfer() = 0;
    };

    static constexpr size_t c_nnull_updates_deactivate = 20;

    /*
     * we should be able to deactivate the load balancing once certain externally
     * specified conditions are met
     */
    size_t m_icycle_active = ~0ul;
    /*
     * total number of blocks across all ranks
     */
    const size_t m_nblock;
    /*
     * number of cycles between reallocation attempts
     */
    const size_t m_period;
    defs::inds m_block_to_rank;
    /*
     * push inward-transferred blocks indices to front
     * pop outward-transferred block indices to back
     */
    std::vector<std::list<size_t>> m_rank_to_blocks;
    std::vector<double> m_mean_work_times;
    std::vector<double> m_gathered_times;
    /*
     * if the least productive rank is working for this proportion of the
     * time worked by the most productive rank, then move a block
     */
    const double m_acceptable_imbalance;
    /*
     * If the load balancing algorithm is stable, a situation should be reached in which
     * the update method finds no reason to reallocate blocks. If this happens for a
     * certain m_nnull_updates_deactivate periods, m_active is automatically set to false.
     */
    size_t m_nnull_updates_deactivate = c_nnull_updates_deactivate;
    size_t m_nnull_updates = 0ul;

private:
    std::list<Dependent*> m_dependents;
    /*
     * for each dependent, append a lambda to a list which is then called in the update method
     */
    std::list<TableBase::recv_cb_t> m_recv_callbacks;

    void refresh_callback_list();

    typename std::list<Dependent*>::iterator add_dependent(Dependent* dependent);

    void erase_dependent(Dependent* dependent);

public:
    RankAllocatorBase(size_t nblock, size_t period, double acceptable_imbalance);

    bool row_mapped_by_dependent(size_t irow){
        for (const auto dep : m_dependents) {
            if (dep->has_row(irow)) return true;
        }
        return false;
    }

    virtual TableBase table() = 0;

    size_t nblock_() const;

    size_t get_nskip_() const;

    void deactivate();

    void activate(size_t icycle);

    bool is_active() const;

    /**
     * checks that the m_block_to_rank map is consistent with its inverse: m_rank_to_blocks
     * @return true if maps are consistent
     */
    bool consistent();

    void update(size_t icycle);

};


template<typename row_t>
class RankAllocator : public RankAllocatorBase {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    MappedTable<row_t>& m_table;
    typedef typename KeyField<row_t>::type key_field_t;

public:
    RankAllocator(MappedTable<row_t>& table, size_t nblock, size_t period, double acceptable_imbalance) :
    RankAllocatorBase(nblock, period, acceptable_imbalance), m_table(table) {}

    TableBase table() override {
        return TableBase(0);
    }

    void record_work_time(const row_t& row, const Timer& work_time) {
        m_mean_work_times[get_block(row)]+=work_time;
    }

    inline size_t get_block(const key_field_t& key) const{
        return key.hash()%m_nblock;
    }

    inline size_t get_block(const row_t& row) const{
        return get_block(KeyField<row_t>::get(row));
    }

    inline size_t get_rank(const key_field_t& key) const{
        return m_block_to_rank[get_block(key)];
    }

    inline size_t get_rank(const row_t& row) const{
        return get_rank(row.m_key_field);
    }
};


#endif //M7_RANKALLOCATOR_H
