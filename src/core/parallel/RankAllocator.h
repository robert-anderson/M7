//
// Created by rja on 27/02/2020.
//

#ifndef M7_RANKALLOCATOR_H
#define M7_RANKALLOCATOR_H

#include "src/core/field/Field.h"
#include "src/core/field/Fields.h"
#include "src/core/table/Table.h"
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

        // called before block transfer is performed
        virtual void before_block_transfer(const defs::inds& irows_send, size_t irank_send, size_t irank_recv) = 0;
        // called after row is inserted into the MappedTable on the recving rank
        virtual void on_row_recv_(size_t irow) = 0;
        // called after block transfer is performed
        virtual void after_block_transfer() = 0;
    };

    static constexpr size_t c_nnull_updates_deactivate = 20;

    Table& m_table;
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
    std::list<Table::recv_cb_t> m_recv_callbacks;

    void refresh_callback_list();

    typename std::list<Dependent*>::iterator add_dependent(Dependent* dependent);

    void erase_dependent(Dependent* dependent);

public:
    RankAllocatorBase(Table& table, size_t nblock, size_t period, double acceptable_imbalance);

    virtual size_t get_block_irow(const size_t& irow) = 0;

    size_t nblock_() const;

    void record_work_time(const size_t& irow, const Timer& work_time);

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


template<typename field_t, typename hash_fn=typename field_t::hash_fn>
class RankAllocator : public RankAllocatorBase {
    typedef typename field_t::view_t view_t;
    field_t& m_field;

public:
    RankAllocator(Table& table, field_t& field, size_t nblock, size_t period, double acceptable_imbalance) :
    RankAllocatorBase(table, nblock, period, acceptable_imbalance), m_field(field)
    {}

private:
    size_t get_block_irow(const size_t &irow) override {
        return get_block(m_field(irow));
    }

public:

    inline size_t get_block(const view_t& key) const{
        return hash_fn()(key)%m_nblock;
    }

    inline size_t get_rank(const view_t& key) const{
        return m_block_to_rank[get_block(key)];
    }
};

namespace ra {
    using Onv = RankAllocator<fields::Onv<>>;
}


#endif //M7_RANKALLOCATOR_H
