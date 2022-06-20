//
// Created by rja on 06/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/basis/Lattice.h"

TEST(Lattice, OrthoObc1d) {
    auto lattice = lattice::make("ortho", {6}, {0});
    std::vector<lattice::adj_row_t> chk_rows = {
            {{1ul, 1}},
            {{0ul, 1}, {2ul, 1}},
            {{1ul, 1}, {3ul, 1}},
            {{2ul, 1}, {4ul, 1}},
            {{3ul, 1}, {5ul, 1}},
            {{4ul, 1}}
    };
    ASSERT_EQ(lattice->m_nsite, chk_rows.size());
    lattice::adj_row_t adj_row;
    for (size_t i=0ul; i<lattice->m_nsite; ++i){
        lattice->get_adj_row(i, adj_row);
        ASSERT_EQ(adj_row, chk_rows[i]);
    }
}

TEST(Lattice, OrthoApbc1d) {
    auto lattice = lattice::make("ortho", {6}, {-1});
    std::vector<lattice::adj_row_t> chk_rows = {
            {{5ul, -1}, {1ul, 1}},
            {{0ul, 1}, {2ul, 1}},
            {{1ul, 1}, {3ul, 1}},
            {{2ul, 1}, {4ul, 1}},
            {{3ul, 1}, {5ul, 1}},
            {{4ul, 1}, {0ul, -1}}
    };
    ASSERT_EQ(lattice->m_nsite, chk_rows.size());
    lattice::adj_row_t adj_row;
    for (size_t i=0ul; i<lattice->m_nsite; ++i){
        lattice->get_adj_row(i, adj_row);
        ASSERT_EQ(adj_row, chk_rows[i]);
    }
}


TEST(Lattice, OrthoObc2d) {
    /*
     *  0  1  2  3
     *  4  5  6  7
     *  8  9  10 11
     */
    auto lattice = lattice::make("ortho", {3, 4}, {0, 0});
    std::vector<lattice::adj_row_t> chk_rows = {
            {{4ul, 1}, {1ul, 1}},
            {{5ul, 1}, {0ul, 1}, {2ul, 1}},
            {{6ul, 1}, {1ul, 1}, {3ul, 1}},
            {{7ul, 1}, {2ul, 1}},
            {{0ul, 1}, {8ul, 1}, {5ul, 1}},
            {{1ul, 1}, {9ul, 1}, {4ul, 1}, {6ul, 1}},
            {{2ul, 1}, {10ul, 1}, {5ul, 1}, {7ul, 1}},
            {{3ul, 1}, {11ul, 1}, {6ul, 1}},
            {{4ul, 1}, {9ul, 1}},
            {{5ul, 1}, {8ul, 1}, {10ul, 1}},
            {{6ul, 1}, {9ul, 1}, {11ul, 1}},
            {{7ul, 1}, {10ul, 1}}
    };
    ASSERT_EQ(lattice->m_nsite, chk_rows.size());
    lattice::adj_row_t adj_row;
    for (size_t i=0ul; i<lattice->m_nsite; ++i){
        lattice->get_adj_row(i, adj_row);
        ASSERT_EQ(adj_row, chk_rows[i]);
    }
}


TEST(Lattice, OrthoPbc2d) {
    /*
     *      +
     *   0  1  2
     * - 3  4  5 -
     *   6  7  8
     *      +
     */
    auto lattice = lattice::make("ortho", {3, 3}, {1, -1});
    std::vector<lattice::adj_row_t> chk_rows = {
            {{6ul, 1}, {3ul, 1}, {2ul, -1}, {1ul, 1}},
            {{7ul, 1}, {4ul, 1}, {0ul, 1}, {2ul, 1}},
            {{8ul, 1}, {5ul, 1}, {1ul, 1}, {0ul, -1}},
            {{0ul, 1}, {6ul, 1}, {5ul, -1}, {4ul, 1}},
            {{1ul, 1}, {7ul, 1}, {3ul, 1}, {5ul, 1}},
            {{2ul, 1}, {8ul, 1}, {4ul, 1}, {3ul, -1}},
            {{3ul, 1}, {0ul, 1}, {8ul, -1}, {7ul, 1}},
            {{4ul, 1}, {1ul, 1}, {6ul, 1}, {8ul, 1}},
            {{5ul, 1}, {2ul, 1}, {7ul, 1}, {6ul, -1}},

    };
    ASSERT_EQ(lattice->m_nsite, chk_rows.size());
    lattice::adj_row_t adj_row;
    for (size_t i=0ul; i<lattice->m_nsite; ++i){
        lattice->get_adj_row(i, adj_row);
        ASSERT_EQ(adj_row, chk_rows[i]);
    }
}