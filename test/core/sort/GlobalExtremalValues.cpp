//
// Created by Robert John Anderson on 2020-08-02.
//

#include "src/core/table/BufferedFields.h"
#include "gtest/gtest.h"
#include "src/core/sort/GlobalExtremalValues.h"

namespace global_extremal_values_test {
    typedef SingleFieldRow<fields::Number<int>> scalar_row_t;
    typedef BufferedTable<scalar_row_t> scalar_table_t;

    static size_t nfind() {
        return 4 * mpi::nrank();
    }

    static size_t get_nrow(size_t irank) {
        return hashing::in_range(irank, 7, 18);
    }

    static std::vector<int> get_data(size_t irank) {
        const auto nrow = get_nrow(irank);
        std::vector<int> out;
        out.reserve(nrow);
        const int vmax = 40;
        for (size_t i = 0ul; i < nrow; ++i)
            out.push_back(-vmax + hashing::in_range((i + 1) * (irank + 1), 0, 2 * vmax));
        return out;
    }

    static void setup(scalar_table_t &table) {
        table.clear();
        auto row = table.m_row;
        const auto data = get_data(mpi::irank());
        for (auto &v: data) {
            row.push_back_jump();
            row.m_field = v;
        }
    }

    static size_t get_global_nrow() {
        size_t nrow = 0ul;
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) nrow += get_nrow(irank);
        return nrow;
    }

    struct Item {
        size_t m_irank, m_irow;
        int m_value;

        Item(size_t irank, size_t irow, int value) : m_irank(irank), m_irow(irow), m_value(value) {}
    };

    static std::vector<Item> get_all_sorted(bool largest, bool absval){
        std::vector<Item> all_items;
        all_items.reserve(get_global_nrow());
        for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
            auto rank_data = get_data(irank);
            for (size_t irow=0ul; irow<rank_data.size(); ++irow)
                all_items.emplace_back(irank, irow, rank_data[irow]);
        }
        std::function<bool(const Item& i1, const Item& i2)> comp_fn;
        if (largest){
            if (absval) comp_fn = [](const Item& i1, const Item& i2){return std::abs(i1.m_value)>std::abs(i2.m_value);};
            else comp_fn = [](const Item& i1, const Item& i2){return i1.m_value>i2.m_value;};
        } else {
            if (absval) comp_fn = [](const Item& i1, const Item& i2){return std::abs(i1.m_value)<std::abs(i2.m_value);};
            else comp_fn = [](const Item& i1, const Item& i2){return i1.m_value<i2.m_value;};
        }

        std::sort(all_items.begin(), all_items.end(), comp_fn);
        return all_items;
    }
}

TEST(GlobalExtremalValues, Nrow) {
    using namespace global_extremal_values_test;
    scalar_table_t table("test", {{}});
    setup(table);
    auto nrow = table.m_hwm;
    nrow = mpi::all_sum(nrow);
    ASSERT_EQ(nrow, get_global_nrow());
}

TEST(GlobalExtremalValues, FindAsc) {
    using namespace global_extremal_values_test;
    scalar_table_t table("test", {{}});
    setup(table);

    for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
        if (mpi::i_am(irank)) {
            std::cout << table.to_string() << std::endl;
        }
        mpi::barrier();
    }

    ASSERT_EQ(1, 2);
    ASSERT_EQ(mpi::all_sum(table.m_hwm), get_global_nrow());

    auto& row = table.m_row;
    GlobalExtremalValues<scalar_row_t, int> gxv(row, row.m_field, false, false, 0ul);
    auto nrow_to_find = nfind();
    gxv.find(nrow_to_find);
    ASSERT_TRUE(mpi::nrank()>0 || gxv.nfound_global()==nrow_to_find);
    ASSERT_GE(gxv.nfound_global(), nrow_to_find);
    auto& sorted_row = gxv.m_global_sorter.m_row;
    if (mpi::i_am_root()) {

        ASSERT_EQ(gxv.m_global_sorter.m_hwm, gxv.nfound_global());
        auto right_order = get_all_sorted(false, false);
        ASSERT_EQ(right_order.size(), get_global_nrow());

        sorted_row.restart();
        for (size_t iitem=0ul; iitem<gxv.nfound_global(); ++iitem){
            const auto& item = right_order[iitem];
            //size_t irank = sorted_row.m_irank;
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            //ASSERT_EQ(irank, item.m_irank);
            sorted_row.step();
        }
        ASSERT_FALSE(sorted_row.in_range());
    }
}


#if 0
TEST(GlobalExtremalValues, WorkingOut) {
    /**
     *
     * suppos
     *
     *
     *
     * given a (asc/desc) sorted vector vi with ni elements for each MPI rank i, select the N (lowest/highest) values
     * over all vi without bulk communication of elements.
     *
     * selection in this case means determining the correct values of indices mi which signify the number of the
     * (lowest/highest) values from rank i to be included in the global set
     *
     * example: we want the lowest 10 values from the following
     * 0  1  9 11 13
     * 1  2  5  6 14 18
     * 1  3  4
     *
     * find divceil(10, 3) = 4 on each rank.
     *
     * expecting to find 12 values
     *
     * we actually only found 11 since rank 2 has only 3 values.
     *
     * does the global (lowest/highest) set exist within the current selection, or do we need to find more local values?
     *
     * the best of the worst (botw) is the (lowest/highest) value among the worst (highest/lowest) values currently found
     * on each rank.
     *
     * in this example, the worst values by rank appear to be [11, 6, 4], but rank 2 has only 3 elements, and so it is
     * not eligible for the botw designation. The point of finding the botw is to find a rank which might have values
     * better than those of the selected set of other ranks.
     *
     * The botw is therefore on rank 1. we can then set the number on each rank that we can be certain belong to the
     * global best. those are the values better (less/greater) than or equal to the botw.
     *
     * These ni are therefore [2, 4, 3] so we have found 9/10 of the desired global values.
     *
     * we then do another iteration to find the missing 1 value(s) there are 2 ranks left active, so find divceil(1, 2)
     * = 1 values on the active ranks 0 and 1.
     *
     * we find additional values [13, 14, -] and can identify rank 0 as having the botw
     *
     * the ni are now [5, 4, 3] so we have found 12/10 of the desired global values.
     *
     * we now
     *
     *
     * vi[mi-1]
     */
    const bool asc = true;
    std::function<bool(const size_t &i1, const size_t &i2)> comp_fn;
    if (asc) comp_fn = [](const size_t &i1, const size_t &i2) { return i1 < i2; };
    else comp_fn = [](const size_t &i1, const size_t &i2) { return i1 > i2; };

    defs::inds ns = {7, 6, 9};
    const size_t vmax = 32;
    std::vector<defs::inds> vs;
    size_t irank = 0ul;
    for (auto n: ns) {
        vs.push_back({});
        for (size_t i = 0ul; i < n; ++i) vs.back().push_back(hashing::in_range(i + n + irank, 0, vmax));
        std::sort(vs.back().begin(), vs.back().end(), comp_fn);
        ++irank;
        std::cout << vs.back() << std::endl;
    }

}

TEST(GlobalExtremalValues, Test) {

    typedef SingleFieldRow <fields::Number<size_t>> row_t;
    BufferedTable <row_t> table("Test", {{}});
    const size_t nrow_per_rank = 40;

    auto get_value = [](size_t irank, size_t ielement) {
        return hashing::in_range((irank + 3) * (ielement + 9), 0, nrow_per_rank * mpi::nrank());
    };

    table.push_back(nrow_per_rank);
    auto row = table.m_row;
    for (row.restart(); row.in_range(); row.step())
        row.m_field = get_value(mpi::irank(), row.index());
    ASSERT_EQ(table.m_hwm, nrow_per_rank);

    const size_t nfind = 10;
    LocalExtremalValues<row_t, size_t, 0ul> lxv(table.m_row, table.m_row.m_field, 1, 1);
    lxv.find(nfind);
    for (size_t i = 0ul; i < lxv.m_xinds.nfound(); ++i) {
        row.jump(lxv.m_xinds[i]);
        std::cout << row.m_field << std::endl;
    }

    /*
    auto row_cmp = row;
    auto cmp_fn = [&](const size_t& irow, const size_t& irow_cmp){
        row.jump(irow);
        row_cmp.jump(irow_cmp);
        return row.m_field < row_cmp.m_field;
    };

    ExtremalIndices xv(cmp_fn);

    const size_t nfind = 10;
    xv.reset(table.m_hwm);
    xv.find(nfind);
    for (size_t i=0ul; i<xv.nfound(); ++i){
        row.jump(xv[i]);
        std::cout << row.m_field << std::endl;
    }
    */
    std::cout << "" << std::endl;
    if (mpi::i_am_root()) {
        defs::inds all_values;
        all_values.reserve(nrow_per_rank * mpi::nrank());
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank)
            for (size_t ielement = 0ul; ielement < nrow_per_rank; ++ielement)
                all_values.emplace_back(get_value(irank, ielement));
        std::sort(all_values.begin(), all_values.end());

        for (size_t i = 0ul; i < nfind; ++i)
            std::cout << all_values[i] << std::endl;
    }

    /*
    const size_t nfind = 8;
    GlobalExtremalValues<field_t> pxv(m_table.m_row.m_field);
    pxv.reset(m_table);

    pxv.find(nfind);
     */

}
#endif