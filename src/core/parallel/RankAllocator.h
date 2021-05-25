/**
 * @file
 * @author Robert John Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 *  FCIQMC and other stochastic wavefunction methods in electronic structure rely on the
 *  ability of the implementation to store large arrays in a distributed manner, with each
 *  rank in the communicator assigned to a subset of the non-zero elements.
 *
 *  Common examples of such large arrays in FCIQMC are the wavefunctions and multidimensional
 *  expectation values such as RDMs and spectral quantities.
 *
 *  The assignment of individual elements to MPI ranks would ultimately defeat the purpose of
 *  the distributed approach, and so elements are assigned their MPI rank block-wise. The block
 *  index is computed via a simple modular hash function acting on the mapped field identifying
 *  the element (ONV for the wavefunction, index array for the MEVs)
 *
 *  Thus, the space of non-zero elements is more or less uniformly spread out across the
 *  compute resources. However their cost (which is of course not a-priori ascertainable) is not.
 *
 *  This cost can be any user-defined measure. In this implementation, the wavefunction uses the
 *  total time elapsed in the spawning loop over occupied ONVs.
 */

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

// forward decl since a ref is required by RankDynamic class
struct RankAllocatorBase;

/**
 * Base class for objects which track rows of a MappedTable whose rows are managed by a
 * RankAllocator
 */
struct RankDynamic {
    RankAllocatorBase &m_ra;
    typename std::list<RankDynamic *>::iterator m_it;

    RankDynamic(RankAllocatorBase &ra);

    ~RankDynamic();

    /**
     * @param irow
     *  row index whose depended-upon status is being determined
     * @return
     *  returns true if the irow-th row in source mapped table is mapped by this RankDynamic
     */
    virtual bool has_row(size_t irow) = 0;

    /**
     * called before physical transfer of rows is performed by MPI
     * @param irows_send
     *  row indices about to be sent
     * @param irank_send
     *  MPI rank on which the rows to be sent are currently stored
     * @param irank_recv
     *  MPI rank to which the rows will be transferred
     */
    virtual void before_block_transfer(const defs::inds &irows_send, size_t irank_send, size_t irank_recv) = 0;

    /**
     * called after each row is inserted into the MappedTable on the recving rank
     * @param irow
     *  row index of last transferred row to be inserted into the MappedTable
     */
    virtual void on_row_recv_(size_t irow) = 0;

    /**
     * called after physical transfer of rows is performed by MPI
     */
    virtual void after_block_transfer() = 0;
};

/**
 * Here we abstract the behavior of the RankAllocator which does not depend on the mapped Field.
 * Objects which reference rows in the mapped table need to be notified if a referenced row is
 * being transferred to another rank, or if they are being given responsibility for tracking a
 * row which is incoming from another rank. Such objects are referred to as RankDynamic, and they
 * are linked to the RankAllocatorBase by a list of addresses.
 */
struct RankAllocatorBase {
    /**
     * The default number of consecutive instances of no change to the ranks to which blocks of
     * elements are allocated ("null updates") before the rank allocation algorithm is turned off.
     */
    static constexpr size_t c_nnull_updates_deactivate = 20;
    /**
     * The "cycle" on which the rank allocator was last made active. This is corresponds to the
     * loop variable of Solver::execute in this program usually, but is not required to.
     */
    size_t m_icycle_active = ~0ul;
    /**
     * total number of blocks across all ranks
     */
    const size_t m_nblock;
    /**
     * number of cycles between reallocation attempts
     */
    const size_t m_period;
    /**
     * A map from block index to MPI rank index
     */
    defs::inds m_block_to_rank;
    /**
     * A map from MPI rank index to a list of all currently allocated blocks
     * push inward-transferred blocks indices to front
     * pop outward-transferred block indices to back
     */
    std::vector<std::list<size_t>> m_rank_to_blocks;
    /**
     * "work time" stands for the figure of merit representing the cost associated with a
     * block. The vector is of length m_nblock but of course only elements at block indices
     * allocated to the local rank are non-zero.
     */
    std::vector<double> m_mean_work_times;
    /**
     * work times must be summed to give the total work time for the local rank, and then
     * that total is gathered in order to determine the "busiest" and "laziest" ranks.
     */
    std::vector<double> m_gathered_total_times;
    /**
     * If the fractional difference in work done by the busiest and laziest ranks is less
     * than this number, then do a null update.
     */
    const double m_acceptable_imbalance;
    /**
     * If the load balancing algorithm is stable, a situation should be reached in which
     * the update method finds no reason to reallocate blocks. If this happens for a
     * certain m_nnull_updates_deactivate periods, m_active is automatically set to false.
     */
    size_t m_nnull_updates_deactivate = c_nnull_updates_deactivate;
    size_t m_nnull_updates = 0ul;

private:
    /**
     * list of addresses of dependent objects which are affected by RankAllocator updates
     */
    std::list<RankDynamic *> m_rank_dynamic_objects;

    /**
     * for each rank-dynamic object, append a lambda to a list which is then called in the update method
     */
    std::list<TableBase::recv_cb_t> m_recv_callbacks;

    /**
     * when rank-dynamic objects are added by the add_dependent method, or removed by the dtor of a
     * rank-dynamic object via the erase_dependent method, the above list of callbacks must be updated
     */
    void refresh_callback_list();

public:
    /**
     * adds a rank-dynamic object to the RankAllocator
     * @param dependent
     *  address of rank-dynamic object to be added to the RankAllocator
     * @return
     *  iterator within list which points to the added rank-dynamic object, this does not change after
     *  calls to erase_dependent for other rank-dynamic objects since the list container is doubly linked.
     */
    typename std::list<RankDynamic *>::iterator add_dependent(RankDynamic *dependent);

    /**
     * removes a rank-dynamic object from the RankAllocator
     * @param dependent
     *  address of rank-dynamic object to be unlinked from this RankAllocator
     */
    void erase_dependent(RankDynamic *dependent);

    RankAllocatorBase(size_t nblock, size_t period, double acceptable_imbalance);

    /**
     * @param irow
     *  index of row whose status as mapped by any added rank-dynamic object is to be determined
     * @return
     *  true if the row is mapped by any rank-dynamic object
     */
    bool row_mapped_by_dependent(size_t irow);

    /**
     * @return
     *  reference to table object in derived class. The RankAllocatorBase makes no reference to
     *  the MappedTable since that would require a template parameter
     */
    virtual TableBase &table() = 0;
    virtual const TableBase &table() const = 0;

    /**
     * @return
     *  number of blocks in the local allocation
     */
    size_t nblock_local() const;

    /**
     * Selective function: only called on the sending rank since the block-wise work
     * times are not shared with all ranks.
     *
     * When selecting a block to transfer, we don't want to send an expensive block
     * as this has the potential to be counterproductive.
     *
     * Thus, send the first block found by iterating from the front of the list with
     * mean_work_time less than the average over all blocks on this MPI rank.
     * @return
     *  number of skips in the list of blocks required to reach the first one which
     *  has a cost less than the average cost over all blocks on this rank.
     */
    size_t get_nskip_() const;

    /**
     * turn off dynamic rank allocation
     */
    void deactivate();

    /**
     * turn on dynamic rank allocation
     * @param icycle
     *  cycle index on which the RankAllocator is activated
     */
    void activate(size_t icycle);

    /**
     * @return
     *  true if dynamic rank allocation is turned on
     */
    bool is_active() const;

    /**
     * checks that the m_block_to_rank map is consistent with its inverse: m_rank_to_blocks
     * @return
     *  true if maps are consistent
     */
    bool consistent();

    /**
     * The update method:
     * 1. decides whether a block transfer should take place
     * 2. if so, determines which rank is busiest and which is laziest
     * 3. via get_nskip_ method, picks a block to transfer from the busiest rank to the
     *    laziest one
     * 4. lets the rank->block and block->rank mappings reflect this reallocation
     * 5. calls the preparatory before_block_transfer virtual method for all dependents
     * 6. calls the block transfer method of the MappedTable to physically communicate the
     *    affected rows
     * 7. calls the after_block_transfer virtual method for all dependents
     * @param icycle
     *  the cycle index, to determine whether or not a block reallocation should be attempted
     */
    void update(size_t icycle);

    virtual size_t get_rank_by_irow(const size_t &irow) const = 0;

    virtual size_t get_block_by_irow(const size_t &irow) const = 0;

};

/**
 * With the business of keeping track of the rank allocation of blocks abstracted to the base
 * class, the templated RankAllocator only needs to couple this functionality to a particular
 * MappedTable
 * @tparam row_t
 *  defined the data layout of the MappedTable, the rank allocation of whose rows is managed by
 *  this class
 */
template<typename row_t>
class RankAllocator : public RankAllocatorBase {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
    MappedTable<row_t> &m_table;
    /**
     * copy the row of m_table to act as a cursor for RankAllocator functions
     */
    row_t m_row;
    /**
     * the block index is computed from a hash function acting on the key field of the row, it
     * must be defined for RankAllocation to work, and compilation will fail here if the key field
     * identifying function is not defined in row_t's class definition (see WalkerTableRow for
     * an example)
     */
    typedef typename KeyField<row_t>::type key_field_t;

public:
    RankAllocator(MappedTable<row_t> &table, size_t nblock, size_t period, double acceptable_imbalance) :
            RankAllocatorBase(nblock, period, acceptable_imbalance), m_table(table), m_row(table.m_row) {}

    /**
     * implements base class method
     * @return
     *  TableBase cast of MappedTable
     */
    const TableBase &table() const override {
        return m_table;
    }

    TableBase &table() override {
        return m_table;
    }

    /**
     * accumulate time spent working on a row to the m_mean_work_times vector
     * @param row
     *  row object pointing to the row whose cost is being recorded
     * @param work_time
     *  time spent on row
     */
    void record_work_time(const row_t &row, const Timer &work_time) {
        m_mean_work_times[get_block(row)] += work_time;
    }

    /**
     * @param key
     *  mapped field instance
     * @return
     *  block index corresponding to mapped field
     */
    inline size_t get_block(const key_field_t &key) const {
        return key.hash() % m_nblock;
    }

    /**
     * @param row
     *  row instance
     * @return
     *  block index corresponding to mapped field of row instance
     */
    inline size_t get_block(const row_t &row) const {
        return get_block(KeyField<row_t>::get(row));
    }

    /**
     * @param key
     *  mapped field instance
     * @return
     *  MPI rank index to which the block corresponding to mapped field is allocated
     */
    inline size_t get_rank(const key_field_t &key) const {
        return m_block_to_rank[get_block(key)];
    }
    /**
     * @param key
     *  row instance
     * @return
     *  MPI rank index to which the block corresponding to mapped field of row object
     *  is allocated
     */
    inline size_t get_rank(const row_t &row) const {
        return get_rank(KeyField<row_t>::get(row));
    }

    /**
     * @param irow
     *  row index to which internal row cursor m_row is moved
     * @return
     *  corresponding MPI rank
     */
    size_t get_rank_by_irow(const size_t &irow) const override {
        m_row.jump(irow);
        return get_rank(m_row);
    }

    /**
     * @param irow
     *  row index to which internal row cursor m_row is moved
     * @return
     *  block index to which row belongs
     */
    size_t get_block_by_irow(const size_t &irow) const override {
        m_row.jump(irow);
        return get_block(m_row);
    }

    /**
     * check that all non-zero key fields in the mapped table are on their allocated rank
     * @return
     *  true if verification passes
     */
    bool verify() const {
        auto row = m_table.m_row;
        const auto& key_field = KeyField<row_t>::get(row);
        for(row.restart(); row.in_range(); row.step()){
            if (key_field.is_zero()) continue;
            if (!mpi::i_am(get_rank(row))) return false;
        }
        return true;
    }

};


#endif //M7_RANKALLOCATOR_H
