//
// Created by rja on 09/07/2021.
//

#ifndef M7_GLOBALEXTREMALVALUES_H
#define M7_GLOBALEXTREMALVALUES_H

#include "LocalExtremalValues.h"
#include "src/core/table/BufferedTable.h"
#include "QuickSorter.h"

template<typename T, size_t nind>
struct GlobalSortingRow : Row {
    fields::Number<size_t> m_irank;
    fields::Numbers<T, nind> m_value;
    GlobalSortingRow(const NdFormat<nind>& format): m_irank(this), m_value(this, format){}
};


template<typename row_t, typename T, size_t nind=0>
struct GlobalExtremalValues {
    typedef LocalExtremalValues<row_t, T, nind> lxv_t;
    typedef GlobalSortingRow<T, nind> global_sort_row_t;
    typedef BufferedTable<global_sort_row_t, false> global_sort_table_t;
    lxv_t m_lxv;
    size_t m_nkeep;
    defs::inds m_global_hwm;
    defs::inds m_global_nfound;
    global_sort_table_t m_global_sorter;
    GlobalExtremalValues(row_t &row, fields::Numbers<T, nind> &field, bool largest, bool absval, size_t ielement_cmp=~0ul):
    m_lxv(row, field, largest, absval, ielement_cmp), m_global_hwm(mpi::nrank()),
    m_global_nfound(mpi::nrank()), m_global_sorter("Global extremal values sorter", {field.m_format}){
        reset();
    }

    size_t get_nrank_with_remainder() const {
        size_t tot = 0ul;
        for (size_t irank=0ul; irank<mpi::nrank(); ++irank) {
            DEBUG_ASSERT_LE(m_global_nfound[irank], m_global_hwm[irank],
                            "number found can't be more than the number of rows stored on the rank!");
            tot+=m_global_nfound[irank]<m_global_hwm[irank];
        }
        return tot;
    }

    void reset() {
        m_lxv.reset();
        mpi::all_gather(m_lxv.hwm(), m_global_hwm);
        m_global_nfound.assign(mpi::nrank(), 0ul);
        m_nkeep = ~0ul;
    }

    size_t nremain() const {
        return m_lxv.nremain();
    }

    size_t nfound_global() const {
        return std::accumulate(m_global_nfound.cbegin(), m_global_nfound.cend(), 0ul);
    }

private:

    /**
     * termination condition for recursive find_required_local_rows
     * @return
     *  true if the eventual set of globally-extremal rows lie within the set of those already found locally
     */
    bool have_all_required_rows() const {
        // TODO
        return true;
//        auto nrank_with_remainder = get_nrank_with_remainder();
//        // if all ranks have found all the rows they store, then there are no more locally-extremal rows to find anywhere
//        if (!nrank_with_remainder) return true;
//        // get the best of the worst "botw"
//        auto irow_local_worst = m_global_nfound[mpi::irank()]-1;
//        auto& row = m_lxv.m_work_row;
//        auto& field = m_lxv.m_work_row_field;
//        row.jump(irow_local_worst);
//        //const T& local_worst = field[m_lxv.m_ielement_cmp];
    }

    void find_required_local_rows(size_t nrow){
        auto nrank_with_remainder = get_nrank_with_remainder();
        if (have_all_required_rows()){
            // no ranks have any more rows to find, or we already have enough
            return;
        }
        DEBUG_ASSERT_LE(nfound_global(), nrow, "should have already returned");
        auto nrow_to_find = nrow - nfound_global();
        size_t nfind_local = std::min(nremain(), integer_utils::divceil(nrow_to_find, nrank_with_remainder));
        m_lxv.find(nfind_local);
        mpi::all_gather(m_lxv.nfound(), m_global_nfound);
        find_required_local_rows(nrow);
    }

    void load_local_values(){
        REQUIRE_TRUE_ALL(nfound_global(), "required local rows haven't been found yet");
        auto& format = m_lxv.m_work_row_field.m_format;
        global_sort_table_t local_loader("Local loader for global extremal values sorter", {format});
        auto nfound_local = m_global_nfound[mpi::irank()];
        static_cast<TableBase&>(local_loader).resize(nfound_local);
        auto& source_row = m_lxv.m_work_row;
        auto& source_field = m_lxv.m_work_row_field;
        auto& loader_row = local_loader.m_row;
        for (size_t i=0ul; i<nfound_local; ++i){
            static_cast<Row&>(source_row).jump(m_lxv[i]);
            static_cast<Row&>(loader_row).push_back_jump();
            loader_row.m_irank = mpi::irank();
            loader_row.m_value = source_field;
        }
        if (mpi::i_am_root())
            static_cast<TableBase&>(m_global_sorter).resize(nfound_global());
        static_cast<TableBase&>(m_global_sorter).gatherv(local_loader);
    }

    void sort(){
        defs::inds nkeep_from_each_rank(mpi::nrank(), 0ul);
        if (mpi::i_am_root()){
            REQUIRE_EQ(m_global_sorter.m_hwm, nfound_global(),
                       "global sorting table should have as many filled rows as total found rows across all ranks");
            auto row1 = m_global_sorter.m_row;
            auto row2 = m_global_sorter.m_row;
            auto cmp_fn = comparators::make_num_field_row_cmp_fn(
                    row1, row1.m_value, row2, row2.m_value, m_lxv.m_value_cmp_fn, m_lxv.m_ielement_cmp);
            QuickSorter qs(cmp_fn);
            qs.reorder_sort(m_global_sorter);
            for (size_t irow=0ul; irow<nfound_global(); ++irow){
                row1.jump(irow);
                ++nkeep_from_each_rank[row1.m_irank];
            }
        }
        mpi::bcast(nkeep_from_each_rank);
        m_nkeep = nkeep_from_each_rank[mpi::irank()];
    }

public:
    void find(size_t nrow){
        find_required_local_rows(nrow);
        load_local_values();
        if (mpi::i_am_root())
            std::cout << m_global_sorter.to_string() << std::endl;
        sort();
        if (mpi::i_am_root())
            std::cout << m_global_sorter.to_string() << std::endl;
    }
};


#endif //M7_GLOBALEXTREMALVALUES_H
