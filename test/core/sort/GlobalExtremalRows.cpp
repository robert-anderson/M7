//
// Created by Robert John Anderson on 2020-08-02.
//

#include "src/core/table/BufferedFields.h"
#include "gtest/gtest.h"
#include "src/core/sort/GlobalExtremalRows.h"

namespace global_extremal_rows_test {
    typedef SingleFieldRow<fields::Number<int>> scalar_row_t;
    typedef BufferedTable<scalar_row_t> scalar_table_t;

    static size_t nfind() {
        return 4 * mpi::nrank();
    }

    static size_t get_nrow(size_t irank) {
        return hashing::in_range(irank, 0, 10);
    }

    static std::vector<int> get_data(size_t irank) {
        const auto nrow = get_nrow(irank);
        std::vector<int> out;
        out.reserve(nrow);
        const int vmax = 20;
        for (size_t i = 0ul; i < nrow; ++i)
            out.push_back(-vmax + hashing::in_range({i, irank}, 0, 2 * vmax));
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

TEST(GlobalExtremalRows, Nrow) {
    using namespace global_extremal_rows_test;
    scalar_table_t table("test", {{}});
    setup(table);
    auto nrow = table.m_hwm;
    nrow = mpi::all_sum(nrow);
    ASSERT_EQ(nrow, get_global_nrow());
}

TEST(GlobalExtremalRows, FindAsc) {
    using namespace global_extremal_rows_test;
    scalar_table_t table("test", {{}});
    setup(table);

    ASSERT_EQ(mpi::all_sum(table.nrow_nonzero()), get_global_nrow());

    auto& row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxr(row, row.m_field, false, false, 0ul);

    auto nrow_to_find = nfind();
    gxr.find(nrow_to_find);

    ASSERT_TRUE(mpi::nrank()>0 || gxr.m_ninclude.m_reduced==nrow_to_find);
    ASSERT_EQ(gxr.m_ninclude.m_reduced, std::min(nrow_to_find, mpi::all_sum(table.nrow_nonzero())));

    /*
     * the global sort only happens on the root rank
     */
    if (mpi::i_am_root()) {
        auto& sorted_row = gxr.m_global_sorter.m_row;
        /*
         * the total number of values considered in the global sorter mus be at least at great as the total number of
         * rows included in the globally extremal set
         */
        ASSERT_GE(gxr.m_global_sorter.m_hwm, gxr.m_ninclude.m_reduced);
        auto right_order = get_all_sorted(false, false);
        ASSERT_EQ(right_order.size(), get_global_nrow());

        sorted_row.restart();
        for (size_t iitem=0ul; iitem<gxr.m_ninclude.m_reduced; ++iitem){
            const auto& item = right_order[iitem];
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            sorted_row.step();
        }
    }
}

