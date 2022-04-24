//
// Created by Oskar Weser on 2022-03-21.
//

#include <gtest/gtest.h>

#include "M7_lib/hamiltonian/frm/SpinSquareFrmHam.h"
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/connection/Connections.h>


TEST(SpinSquareFrmHam, Elements) {
    SpinSquareFrmHam s2({8, 6, 0});
    {
        buffered::FrmOnv det(s2.m_hs);
        det = {{0, 1, 2, 3}, {0, 1, 2, 4}};

        ASSERT_EQ(s2.get_element_0000(det), 3);

        conn::FrmOnv conn(det);
        conn.m_ann.set({0, 3});
        conn.m_cre.set({0, 5});
        ASSERT_EQ(s2.get_element_1100(det, conn), 0.0);

        conn.m_ann.set({0, 3}, {1, 4});
        conn.m_cre.set({0, 4}, {1, 3});
        ASSERT_EQ(s2.get_element_2200(det, conn), 1.0);

        conn.m_ann.set({0, 2}, {0, 3});
        conn.m_cre.set({0, 4}, {0, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);

        conn.m_ann.set({0, 1}, {1, 2});
        conn.m_cre.set({0, 5}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);

        conn.m_ann.set({0, 1}, {1, 1});
        conn.m_cre.set({0, 5}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);

        conn.m_ann.set({0, 3}, {1, 4});
        conn.m_cre.set({0, 5}, {1, 3});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);

        conn.m_ann.set({1, 1}, {1, 2});
        conn.m_cre.set({1, 3}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);
    }

    {
        buffered::FrmOnv det(s2.m_hs);
        det = {{0, 1, 2, 4}, {0, 1, 2, 3}};

        ASSERT_EQ(s2.get_element_0000(det), 3);

        conn::FrmOnv conn(det);
        conn.m_ann.set({0, 4});
        conn.m_cre.set({0, 5});
        ASSERT_EQ(s2.get_element_1100(det, conn), 0.0);

        conn.m_ann.set({0, 4}, {1, 3});
        conn.m_cre.set({0, 3}, {1, 4});
        ASSERT_EQ(s2.get_element_2200(det, conn), 1.0);

        conn.m_ann.set({0, 4}, {1, 3});
        conn.m_cre.set({0, 3}, {1, 4});
        ASSERT_EQ(s2.get_element_2200(det, conn), 1.0);

        conn.m_ann.set({0, 2}, {0, 4});
        conn.m_cre.set({0, 3}, {0, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);

        conn.m_ann.set({1, 1}, {1, 2});
        conn.m_cre.set({1, 4}, {1, 5});
        ASSERT_EQ(s2.get_element_2200(det, conn), 0.0);
    }
}
