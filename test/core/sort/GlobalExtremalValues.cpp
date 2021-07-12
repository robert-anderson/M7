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

    for (size_t irank=0ul; irank<mpi::nrank(); ++irank){
        if (mpi::i_am(irank)) {
            std::cout << table.to_string() << std::endl;
        }
        mpi::barrier();
    }

    ASSERT_EQ(1, 2);
    ASSERT_EQ(mpi::all_sum(table.m_hwm), get_global_nrow());

    auto& row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxv(row, row.m_field, false, false, 0ul);
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
