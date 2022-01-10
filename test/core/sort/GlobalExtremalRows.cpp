//
// Created by Robert John Anderson on 2020-08-02.
//

#include "src/core/table/BufferedFields.h"
#include "gtest/gtest.h"
#include "src/core/sort/GlobalExtremalRows.h"

namespace global_extremal_rows_test {
    typedef SingleFieldRow<field::Number<int>> scalar_row_t;
    typedef BufferedTable<scalar_row_t> scalar_table_t;

    static size_t nfind() {
        return 6 * mpi::nrank();
    }

    static size_t get_nrow(size_t irank) {
        return hashing::in_range(irank, 0, 26);
    }

    static bool skip(size_t irank, size_t i){
        // clear one row in every 10
        return !hashing::in_range({irank, i}, 0, 10);
    }

    static size_t get_nrow_nonzero(size_t irank) {
        size_t nrow = 0ul;
        for (size_t irow=0ul; irow<get_nrow(irank); ++irow)
            if (!skip(irank, irow)) ++nrow;
        return nrow;
    }

    static std::vector<int> get_data(size_t irank) {
        const auto nrow = get_nrow(irank);
        std::vector<int> out;
        out.reserve(nrow);
        const int vmax = 20;
        for (size_t i = 0ul; i < nrow; ++i)
            out.push_back(-vmax + int(hashing::in_range({i, irank}, 0, 2 * vmax)));
        return out;
    }

    static std::vector<int> get_data_with_skips(size_t irank) {
        auto data = get_data(irank);
        auto out = data;
        out.clear();
        for (size_t i =0ul; i<data.size(); ++i) if (!skip(irank, i)) out.push_back(data[i]);
        return out;
    }

    static void setup(scalar_table_t &table) {
        table.clear();
        auto row = table.m_row;
        const auto data = get_data(mpi::irank());
        table.resize(data.size());
        for (auto &v: data) {
            row.push_back_jump();
            if (skip(mpi::irank(), row.index())) table.clear(row.index());
            else row.m_field = v;
        }
    }

    static size_t get_global_nrow() {
        size_t nrow = 0ul;
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) nrow += get_nrow(irank);
        return nrow;
    }

    static size_t get_global_nrow_nonzero() {
        size_t nrow = 0ul;
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) nrow += get_nrow_nonzero(irank);
        return nrow;
    }

    struct Item {
        size_t m_irank, m_irow;
        int m_value;

        Item(size_t irank, size_t irow, int value) : m_irank(irank), m_irow(irow), m_value(value) {}
    };

    static std::vector<Item> get_all_sorted(bool largest, bool absval) {
        std::vector<Item> all_items;
        all_items.reserve(get_global_nrow());
        for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            auto rank_data = get_data(irank);
            for (size_t irow = 0ul; irow < rank_data.size(); ++irow){
                if (!skip(irank, irow)) all_items.emplace_back(irank, irow, rank_data[irow]);
            }
        }
        std::function<bool(const Item &i1, const Item &i2)> comp_fn;
        if (largest) {
            if (absval)
                comp_fn = [](const Item &i1, const Item &i2) {
                    return std::abs(i1.m_value) > std::abs(i2.m_value);
                };
            else comp_fn = [](const Item &i1, const Item &i2) { return i1.m_value > i2.m_value; };
        } else {
            if (absval)
                comp_fn = [](const Item &i1, const Item &i2) {
                    return std::abs(i1.m_value) < std::abs(i2.m_value);
                };
            else comp_fn = [](const Item &i1, const Item &i2) { return i1.m_value < i2.m_value; };
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

    auto &row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxr(row, row.m_field, false, false, 0ul);

    auto nrow_to_find = nfind();
    gxr.find(nrow_to_find);

    ASSERT_TRUE(mpi::nrank() > 0 || gxr.m_ninclude.m_reduced == nrow_to_find);
    ASSERT_EQ(gxr.m_ninclude.m_reduced, std::min(nrow_to_find, mpi::all_sum(table.nrow_nonzero())));

    /*
     * the global sort only happens on the root rank
     */
    if (mpi::i_am_root()) {
        auto &sorted_row = gxr.m_global_sorter.m_row;
        /*
         * the total number of values considered in the global sorter mus be at least at great as the total number of
         * rows included in the globally extremal set
         */
        ASSERT_GE(gxr.m_global_sorter.m_hwm, gxr.m_ninclude.m_reduced);
        auto right_order = get_all_sorted(false, false);
        ASSERT_EQ(right_order.size(), get_global_nrow_nonzero());

        sorted_row.restart();
        for (size_t iitem = 0ul; iitem < gxr.m_ninclude.m_reduced; ++iitem) {
            const auto &item = right_order[iitem];
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            sorted_row.step();
        }
    }
}


TEST(GlobalExtremalRows, FindDesc) {
    using namespace global_extremal_rows_test;
    scalar_table_t table("test", {{}});
    setup(table);

    auto &row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxr(row, row.m_field, true, false, 0ul);

    auto nrow_to_find = nfind();
    gxr.find(nrow_to_find);

    ASSERT_TRUE(mpi::nrank() > 0 || gxr.m_ninclude.m_reduced == nrow_to_find);
    ASSERT_EQ(gxr.m_ninclude.m_reduced, std::min(nrow_to_find, mpi::all_sum(table.nrow_nonzero())));

    /*
     * the global sort only happens on the root rank
     */
    if (mpi::i_am_root()) {
        auto &sorted_row = gxr.m_global_sorter.m_row;
        /*
         * the total number of values considered in the global sorter mus be at least at great as the total number of
         * rows included in the globally extremal set
         */

        ASSERT_GE(gxr.m_global_sorter.m_hwm, gxr.m_ninclude.m_reduced);
        auto right_order = get_all_sorted(true, false);
        ASSERT_EQ(right_order.size(), get_global_nrow_nonzero());

        sorted_row.restart();
        for (size_t iitem = 0ul; iitem < gxr.m_ninclude.m_reduced; ++iitem) {
            const auto &item = right_order[iitem];
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            sorted_row.step();
        }
    }
}