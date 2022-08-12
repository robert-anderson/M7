//
// Created by rja on 06/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/basis/Lattice.h"

namespace lattice_test {
    typedef sparse::MatrixElement<int> elem_t;
    typedef v_t<elem_t> row_t;
    typedef v_t<row_t> rows_t;

    row_t make_row(row_t::const_iterator cbegin, row_t::const_iterator cend) {
        return {cbegin, cend};
    }
    bool verify_rows(const lattice::SubLattice& lattice, const rows_t& rows){
        const auto& adj = lattice.m_sparse_adj;
        if (lattice.m_nsite != rows.size()) return false;
        for (uint_t isite=0ul; isite < lattice.m_nsite; ++isite){
            auto row = make_row(adj.cbegin(isite), adj.cend(isite));
            if (row != rows[isite]) return false;
        }
        return true;
    }
}

TEST(SubLattice, OrthoObc1d) {
    using namespace lattice_test;
    auto lattice = lattice::make("ortho", {6}, {0});
    rows_t chk_rows = {
        {{1ul, 1}},
        {{0ul, 1}, {2ul, 1}},
        {{1ul, 1}, {3ul, 1}},
        {{2ul, 1}, {4ul, 1}},
        {{3ul, 1}, {5ul, 1}},
        {{4ul, 1}}
    };
    ASSERT_TRUE(verify_rows(*lattice, chk_rows));
}

TEST(SubLattice, OrthoApbc1d) {
    using namespace lattice_test;
    auto lattice = lattice::make("ortho", {6}, {-1});
    rows_t chk_rows = {
            {{5ul, -1}, {1ul, 1}},
            {{0ul, 1}, {2ul, 1}},
            {{1ul, 1}, {3ul, 1}},
            {{2ul, 1}, {4ul, 1}},
            {{3ul, 1}, {5ul, 1}},
            {{4ul, 1}, {0ul, -1}}
    };
    ASSERT_TRUE(verify_rows(*lattice, chk_rows));
}


TEST(SubLattice, OrthoObc2d) {
    /*
     *  0  1  2  3
     *  4  5  6  7
     *  8  9  10 11
     */
    using namespace lattice_test;
    auto lattice = lattice::make("ortho", {3, 4}, {0, 0});
    rows_t chk_rows = {
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
    ASSERT_TRUE(verify_rows(*lattice, chk_rows));
}

TEST(Lattice, OrthoPbc2d) {
    /*
     *      +
     *   0  1  2
     * - 3  4  5 -
     *   6  7  8
     *      +
     */
    using namespace lattice_test;
    auto lattice = lattice::make("ortho", {3, 3}, {1, -1});
    rows_t chk_rows = {
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
    ASSERT_TRUE(verify_rows(*lattice, chk_rows));
}


TEST(Lattice, NullNextNearest) {
    using namespace lattice_test;
    auto lattice = lattice::make();
    rows_t chk_rows;
    const auto nn = lattice->make_next_nearest();
    ASSERT_TRUE(verify_rows(*nn, chk_rows));
}

TEST(Lattice, OrthoObc1dNextNearest) {
    using namespace lattice_test;
    auto lattice = lattice::make("ortho", {6}, {0});
    rows_t chk_rows = {
            {{2ul, 1}},
            {{3ul, 1}},
            {{0ul, 1}, {4ul, 1}},
            {{1ul, 1}, {5ul, 1}},
            {{2ul, 1}},
            {{3ul, 1}}
    };
    const auto nn = lattice->make_next_nearest();
    ASSERT_TRUE(verify_rows(*nn, chk_rows));
}
