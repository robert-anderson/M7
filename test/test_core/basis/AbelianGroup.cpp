//
// Created by Robert J. Anderson on 07/06/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/basis/AbelianGroup.h"

TEST(AbelianGroup, D2h){
    strv_t labels = {"A1g", "B1g", "B2g", "B3g", "A1u", "B1u", "B2u", "B3u"};
    AbelianGroup d2h(labels, [](const uint_t& iirrep, const uint_t& jirrep){return iirrep^jirrep;});
    /*
     * D2h A1g B1g B2g B3g A1u B1u B2u B3u
     * A1g A1g
     * B1g B1g A1g
     * B2g B2g B3g A1g
     * B3g B3g B2g B1g A1g
     * A1u A1u B1u B2u B3u A1g
     * B1u B1u A1u B3u B2u B1g A1g
     * B2u B2u B3u A1u B1u B2g B3g A1g
     * B3u B3u B2u B1u A1u B3g B2g B1g A1g
     */
    enum d2h_inds {A1g, B1g, B2g, B3g, A1u, B1u, B2u, B3u};
    for (uint_t iirrep = A1g; iirrep<=B3u; iirrep++){
        ASSERT_EQ(d2h.product(iirrep, iirrep), A1g);
        ASSERT_EQ(d2h.product(iirrep, A1g), iirrep);
    }
    uint_t iirrep = B1g;
    ASSERT_EQ(d2h.product(iirrep, B2g), B3g);
    ASSERT_EQ(d2h.product(iirrep, B3g), B2g);
    ASSERT_EQ(d2h.product(iirrep, A1u), B1u);
    ASSERT_EQ(d2h.product(iirrep, B1u), A1u);
    ASSERT_EQ(d2h.product(iirrep, B2u), B3u);
    ASSERT_EQ(d2h.product(iirrep, B3u), B2u);
}