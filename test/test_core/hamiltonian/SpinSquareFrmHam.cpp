//
// Created by Oskar Weser on 2022-03-21.
//

#include <gtest/gtest.h>

#include <M7_lib/hamiltonian/SpinSquareFrmHam.h>
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/connection/Connections.h>


TEST(SpinSquareFrmHam, TestSpinHamiltonian) {
    SpinSquareFrmHam s2(8, 6, 2);

    {
        buffered::FrmOnv det_i(s2.m_nsite);
        det_i = {{0, 1, 2, 3}, {0, 1, 2, 4}};

        ASSERT_EQ(s2.get_element_0000(det_i), 3);

        conn::FrmOnv conn(det_i);
        conn.set({0, 3}, {0, 5});
        ASSERT_EQ(s2.get_element_1100(det_i, conn), 0.0);

        conn.set({0, 3}, {1, 4}, {0, 4}, {1, 3});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 1.0);

        conn.set({0, 2}, {0, 3}, {0, 4}, {0, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);

        conn.set({0, 1}, {1, 2}, {0, 5}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);

        conn.set({0, 1}, {1, 1}, {0, 5}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);

        conn.set({0, 3}, {1, 4}, {0, 5}, {1, 3});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);

        conn.set({1, 1}, {1, 2}, {1, 3}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);
    }

    {
        buffered::FrmOnv det_i(s2.m_nsite);
        det_i = {{0, 1, 2, 4}, {0, 1, 2, 3}};

        ASSERT_EQ(s2.get_element_0000(det_i), 3);

        conn::FrmOnv conn(det_i);
        conn.set({0, 4}, {0, 5});
        ASSERT_EQ(s2.get_element_1100(det_i, conn), 0.0);

        conn.set({0, 4}, {1, 3}, {0, 3}, {1, 4});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 1.0);

        conn.set({0, 4}, {1, 3}, {0, 3}, {1, 4});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 1.0);

        conn.set({0, 2}, {0, 4}, {0, 3}, {0, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);

        conn.set({1, 1}, {1, 2}, {1, 4}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det_i, conn), 0.0);
    }
}
