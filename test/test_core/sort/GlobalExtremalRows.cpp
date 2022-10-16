//
// Created by Robert John Anderson on 2020-08-02.
//

#include "M7_lib/table/BufferedFields.h"
#include "gtest/gtest.h"
#include "M7_lib/sort/GlobalExtremalRows.h"

namespace global_extremal_rows_test {
    typedef SingleFieldRow<field::Number<int>> scalar_row_t;
    typedef buffered::Table<scalar_row_t> scalar_table_t;

    static uint_t nfind() {
        return 6 * mpi::nrank();
    }

    static uint_t get_nrow(uint_t irank) {
        return hash::in_range(irank, 0, 26);
    }

    static bool skip(uint_t irank, uint_t i){
        // clear one row in every 10
        return !hash::in_range({irank, i}, 0, 10);
    }

    static uint_t get_nrow_nonzero(uint_t irank) {
        uint_t nrow = 0ul;
        for (uint_t irow=0ul; irow<get_nrow(irank); ++irow)
            if (!skip(irank, irow)) ++nrow;
        return nrow;
    }

    static v_t<int> get_data(uint_t irank) {
        const auto nrow = get_nrow(irank);
        v_t<int> out;
        out.reserve(nrow);
        const int vmax = 20;
        for (uint_t i = 0ul; i < nrow; ++i)
            out.push_back(-vmax + int(hash::in_range({i, irank}, 0, 2 * vmax)));
        return out;
    }

    static v_t<int> get_data_with_skips(uint_t irank) {
        auto data = get_data(irank);
        auto out = data;
        out.clear();
        for (uint_t i =0ul; i<data.size(); ++i) if (!skip(irank, i)) out.push_back(data[i]);
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

    static uint_t get_global_nrow() {
        uint_t nrow = 0ul;
        for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) nrow += get_nrow(irank);
        return nrow;
    }

    static uint_t get_global_nrow_nonzero() {
        uint_t nrow = 0ul;
        for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) nrow += get_nrow_nonzero(irank);
        return nrow;
    }

    struct Item {
        uint_t m_irank, m_irow;
        int m_value;

        Item(uint_t irank, uint_t irow, int value) : m_irank(irank), m_irow(irow), m_value(value) {}
    };

    static v_t<Item> get_all_sorted(bool largest, bool absval) {
        v_t<Item> all_items;
        all_items.reserve(get_global_nrow());
        for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
            auto rank_data = get_data(irank);
            for (uint_t irow = 0ul; irow < rank_data.size(); ++irow){
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
    scalar_table_t table("test", {});
    setup(table);
    auto nrow = table.nrow_in_use();
    nrow = mpi::all_sum(nrow);
    ASSERT_EQ(nrow, get_global_nrow());
}

TEST(GlobalExtremalRows, FindAsc) {
    using namespace global_extremal_rows_test;
    scalar_table_t table("test", {});
    setup(table);

    auto &row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxr(row, row.m_field, false, false, 0ul);

    auto nrow_to_find = nfind();
    gxr.find(nrow_to_find);

    ASSERT_TRUE(mpi::nrank() > 0 || gxr.m_ninclude.m_reduced == nrow_to_find);
    ASSERT_EQ(gxr.m_ninclude.m_reduced, std::min(nrow_to_find, mpi::all_sum(table.nrecord())));

    /*
     * the global sort only happens on the root rank
     */
    if (mpi::i_am_root()) {
        auto &sorted_row = gxr.m_global_sorter.m_row;
        /*
         * the total number of values considered in the global sorter mus be at least at great as the total number of
         * rows included in the globally extremal set
         */
        ASSERT_GE(gxr.m_global_sorter.nrow_in_use(), gxr.m_ninclude.m_reduced);
        auto right_order = get_all_sorted(false, false);
        ASSERT_EQ(right_order.size(), get_global_nrow_nonzero());

        sorted_row.restart();
        for (uint_t iitem = 0ul; iitem < gxr.m_ninclude.m_reduced; ++iitem) {
            const auto &item = right_order[iitem];
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            sorted_row.step();
        }
    }
}


TEST(GlobalExtremalRows, FindDesc) {
    using namespace global_extremal_rows_test;
    scalar_table_t table("test", {});
    setup(table);

    auto &row = table.m_row;
    GlobalExtremalRows<scalar_row_t, int> gxr(row, row.m_field, true, false, 0ul);

    auto nrow_to_find = nfind();
    gxr.find(nrow_to_find);

    ASSERT_TRUE(mpi::nrank() > 0 || gxr.m_ninclude.m_reduced == nrow_to_find);
    ASSERT_EQ(gxr.m_ninclude.m_reduced, std::min(nrow_to_find, mpi::all_sum(table.nrecord())));

    /*
     * the global sort only happens on the root rank
     */
    if (mpi::i_am_root()) {
        auto &sorted_row = gxr.m_global_sorter.m_row;
        /*
         * the total number of values considered in the global sorter mus be at least at great as the total number of
         * rows included in the globally extremal set
         */

        ASSERT_GE(gxr.m_global_sorter.nrow_in_use(), gxr.m_ninclude.m_reduced);
        auto right_order = get_all_sorted(true, false);
        ASSERT_EQ(right_order.size(), get_global_nrow_nonzero());

        sorted_row.restart();
        for (uint_t iitem = 0ul; iitem < gxr.m_ninclude.m_reduced; ++iitem) {
            const auto &item = right_order[iitem];
            int value = sorted_row.m_value;
            ASSERT_EQ(value, item.m_value);
            sorted_row.step();
        }
    }
}