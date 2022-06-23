//
// Created by Robert J. Anderson on 09/07/2021.
//

#ifndef M7_GLOBALEXTREMALROWS_H
#define M7_GLOBALEXTREMALROWS_H

#include <M7_lib/parallel/Reduction.h>
#include <M7_lib/table/BufferedTable.h>

#include "LocalExtremalRows.h"
#include "LambdaQuickSorter.h"

/**
 * A row for the loading and gathering of locally sorted values and their MPI rank index of origin
 * @tparam T
 *  type of value to be sorted on
 */
template<typename T>
struct GlobalSortingRow : Row {
    field::Number<size_t> m_irank;
    field::Number<T> m_value;

    GlobalSortingRow() : m_irank(this), m_value(this) {}
};

/**
 * applies the LocalExtremalRows class to
 * @tparam row_t
 * @tparam T
 * @tparam nind
 */
template<typename row_t, typename T, size_t nind = 0>
struct GlobalExtremalRows {
    typedef LocalExtremalRows<row_t, T, nind> lxr_t;
    typedef GlobalSortingRow<T> global_sort_row_t;
    typedef BufferedTable<global_sort_row_t, false> global_sort_table_t;
    /**
     * LocalExtremalRows for determining the best values on each process independently
     */
    lxr_t m_lxr;
    /**
     * the find method will typically require the finding of a number of LXRs greater than the number of best values
     * on this rank that are included in the GXR set. the local value of the reducible is always <= the return value of
     * m_lxv.nfound(), and if it is non-zero, the LXR at index m_local-1 is the "worst of the best" values on the local
     * MPI rank.
     *
     * After the sum reduction, the m_reduced member will hold the total number of included rows over all ranks
     *
     * After the global sort, the local number is revised lower or kept the same.
     */
    Reduction<size_t> m_ninclude;
    /**
     * once the number of included rows over all ranks is at least the number requested, the rank index and sorting
     * value will be loaded into a local table, and then gathered into this table on the root rank only, it is then
     * sorted and the refined values of m_ninclude are scattered back to all ranks in the MPI communicator
     */
    global_sort_table_t m_global_sorter;

    GlobalExtremalRows(row_t &row, field::Numbers<T, nind> &field, bool largest, bool absval, defs::ivec_t inds_to_cmp) :
            m_lxr(row, field, largest, absval, inds_to_cmp),
            m_global_sorter("Global extremal rows sorter", {{}}) {
        reset();
    }

    GlobalExtremalRows(row_t &row, field::Numbers<T, nind> &field, bool largest, bool absval, size_t ind_to_cmp) :
            GlobalExtremalRows(row, field, largest, absval, defs::ivec_t{ind_to_cmp}){}

    /**
     * @return
     *  the number of MPI ranks with any nonzero rows left to find in partial sorting operations
     */
    size_t get_nrank_with_remainder() const {
        size_t count = m_lxr.nremain() > 0;
        return mpi::all_sum(count);
    }
    /**
     * start the sort from scratch
     */
    void reset() {
        m_lxr.reset();
        m_ninclude.m_local = 0ul;
    }

private:
    /**
     * in response to the newly expanded locally-extremal sets of rows, determine which of these should be included in
     * the global set by identifying the rank with the "best of the worst values", If this rank still has nonzero rows
     * remaining, its worst value is the criterion for inclusion of rows from all ranks.
     */
    void update_ninclude() {
        /*
         * get the local worst value
         */
        T local_worst_value = 0;
        if (m_lxr.nfound()) local_worst_value = m_lxr.get_value(m_lxr.nfound()-1);
        std::vector<T> local_worst_values(mpi::nrank());
        std::vector<size_t> local_nfounds(mpi::nrank());
        std::vector<size_t> local_nrow_nonzero(mpi::nrank());
        mpi::all_gather(local_worst_value, local_worst_values);
        mpi::all_gather(m_lxr.nfound(), local_nfounds);
        mpi::all_gather(m_lxr.nrow_nonzero(), local_nrow_nonzero);
        defs::ivec_t ignored_ranks;
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            if (local_nfounds[irank]==local_nrow_nonzero[irank]) ignored_ranks.push_back(irank);
        }
        /*
         * set the number included to the maximum value and then decrease till the worst included value on each rank is
         * as good as or better than the botw values
         */
        m_ninclude.m_local = m_lxr.nfound();
        if (ignored_ranks.size() < mpi::nrank()) {
            /*
             * setup a function to order the worst sorted values from each rank
             */
            auto index_cmp_fn = [&](const size_t &irank, const size_t &jrank) {
                return m_lxr.m_value_cmp_fn(local_worst_values[irank], local_worst_values[jrank]);
            };
            ExtremalIndices xi(index_cmp_fn);
            xi.reset(mpi::nrank(), ignored_ranks);
            // find the "best of the worst"
            xi.find(1);
            if (xi.nfound()) {
                auto irank_botw = xi[0];
                auto botw = local_worst_values[irank_botw];
                while (m_ninclude.m_local && m_lxr.cmp_values(botw, m_lxr.get_value(m_ninclude.m_local-1))) {
                    /*
                     * reduce the number of included rows until:
                     *  there are none left, or
                     *  the best of the worst values is no longer better than the worst value on this rank
                     */
                    --m_ninclude.m_local;
                }
            }
        }
        m_ninclude.all_sum();
    }
    /**
     * find more locally extremal rows. this is done so as to get close to the target number of rows but not overshoot
     * by too much, since we eventually need to communicate all included values to the root rank for sorting
     * @param nrow
     *  target number of row to find globally
     */
    void find_required_local_rows(size_t nrow) {
        REQUIRE_NE_ALL(m_ninclude.m_reduced, ~0ul, "required local row reduction hasn't been performed");
        auto nrank_with_remainder = get_nrank_with_remainder();
        if (nrank_with_remainder == 0 || m_ninclude.m_reduced >= nrow) {
            // no ranks have any more rows to find or we already have enough rows
            return;
        }
        /*
         * the general situation on each rank can be demonstrated as follows
         * I: included rows from the locally extreme rows found so far
         * X: the rest of the locally extreme rows (excluded)
         * O: the rest of the nonzero rows on the rank which are not in any particular order
         * IIIIIXXXXXXXXXOOOOOOOOOOO|
         * e.g. a call to m_lxv.find(4) would change this picture to
         * IIIIIXXXXXXXXXXXXXOOOOOOO|
         * the number of rows we want to make available for comparison in the next iteration is:
         */
        auto nfind_local = integer::divceil(nrow - m_ninclude.m_reduced, nrank_with_remainder);
        /*
         * but this would find new locally extreme rows from X/O threshold, whereas we really want to find these from
         * the I/X threshold, since it doesn't make sense to find more rows on a process with few already included, thus
         * the offset must be subtracted. It may be that there are enough X rows on the rank that no local find is even
         * needed
         */
        nfind_local -= std::min(nfind_local, m_lxr.nfound() - m_ninclude.m_local);
        /*
         * now find the required local extremal rows
         */
        m_lxr.find(nfind_local);
        /*
         * and update the best of the worst "I" rows
         */
        update_ninclude();
        /*
         * call recursively
         */
        find_required_local_rows(nrow);
    }
    /**
     * once the required number of locally extremal rows has been found, the next step is to prepare for the sorting
     * operation by loading the locally extremal included values along with their rank of origin into a "local_loader"
     * table, whereupon these tables are gathered onto the root rank into a table with the same row type
     *
     * this could of course have been done with lower communication overhead by just sending the values and associating
     * the offset of the gathered data with the rank index. but crucially we need to keep track of the ranks of origin
     * of each value in the sort so we know how many to include from each rank in the final globally extremal set, and
     * so the rank-value rowed table would need to be constructed for the global sort anyway. The implemented solution
     * is deemed to be the cleanest at this moment, but perhaps the communication overhead could be reduced as described
     * here in future if it proves to be problematic.
     */
    void load_values_for_sorting() {
        REQUIRE_NE_ALL(m_ninclude.m_reduced, ~0ul, "required local row reduction hasn't been performed");
        REQUIRE_TRUE_ALL(m_ninclude.m_reduced, "required local rows haven't been found yet");
        global_sort_table_t local_loader("Local loader for global extremal rows sorter", {{}});
        static_cast<TableBase &>(local_loader).resize(m_ninclude.m_local);
        auto &source_row = m_lxr.m_work_row;
        auto &source_field = m_lxr.m_work_row_field;
        auto &loader_row = local_loader.m_row;
        for (size_t i = 0ul; i < m_ninclude.m_local; ++i) {
            static_cast<Row &>(source_row).jump(m_lxr[i]);
            static_cast<Row &>(loader_row).push_back_jump();
            loader_row.m_irank = mpi::irank();
            loader_row.m_value = source_field.sum_over(m_lxr.m_inds_to_cmp);
        }
        if (mpi::i_am_root())
            static_cast<TableBase &>(m_global_sorter).resize(m_ninclude.m_reduced);
        static_cast<TableBase &>(m_global_sorter).gatherv(local_loader);
    }
    /**
     * do an inplace quicksort on the gathered globally extremal rank-value pairs and loop over the result keeping tally
     * of the number of times each rank of origin is encountered before the number of rows to be included is reached
     * @param nrow
     *  desired number of rows to find
     */
    void sort(size_t nrow) {
        nrow = std::min(nrow, m_ninclude.m_reduced);
        defs::ivec_t ninclude_each_rank(mpi::nrank(), 0ul);
        if (mpi::i_am_root()) {
            REQUIRE_EQ(m_global_sorter.m_hwm, m_ninclude.m_reduced,
                       "global sorting table should have as many filled rows as total found rows across all ranks");
            auto row1 = m_global_sorter.m_row;
            auto row2 = m_global_sorter.m_row;
            /*
             * values are averaged into one in the global_sort_table_t tables, so just one index to compare;
             */
            auto cmp_fn = comparators::make_num_field_row_cmp_fn(
                    row1, row1.m_value, row2, row2.m_value, m_lxr.m_value_cmp_fn, {0ul});
            LambdaQuickSorter qs(cmp_fn);
            qs.reorder_sort(m_global_sorter);
            for (size_t irow = 0ul; irow < nrow; ++irow) {
                row1.jump(irow);
                ++ninclude_each_rank[row1.m_irank];
            }
        }
        mpi::bcast(ninclude_each_rank);
        m_ninclude.m_local = ninclude_each_rank[mpi::irank()];
        m_ninclude.all_sum();
    }

public:
    /**
     * exposes a public interface for the private methods which compute the globally-extreme rows
     * @param nrow
     *  number of rows to find across all ranks, this will
     */
    void find(size_t nrow) {
        update_ninclude();
        if (!mpi::all_sum(m_lxr.m_table.nrow_nonzero())) return;
        find_required_local_rows(nrow);
        load_values_for_sorting();
        sort(nrow);
    }
    /**
     * @param i
     *  ordinal index of the row among the local inclusions in the globally extremal set. with i=0, the row index
     *  returned is the best value on this rank, and with i=m_ninclude-1 (provided that m_ninclude>0) the returned
     *  index is that of the worst value on this rank to be included in the globally extremal set
     * @return
     *  the row index associated with the ordinal index
     */
    const size_t& operator[](const size_t& i) const {
        DEBUG_ASSERT_LT(i, m_ninclude.m_local, "the specified index was not included in the globally extremal set");
        return m_lxr[i];
    }

    void gatherv(Table<row_t>& dst, size_t iroot=0ul) const {
        BufferedTable<row_t> m_local("locally included globally extreme rows", m_lxr.m_work_row);
        m_local.push_back(m_ninclude.m_local);
        for (size_t iinclude=0ul; iinclude<m_ninclude.m_local; ++iinclude){
            m_local.copy_row_in(m_lxr.m_table, (*this)[iinclude], iinclude);
        }
        dst.gatherv(m_local, iroot);
    }
};


#endif //M7_GLOBALEXTREMALROWS_H
