//
// Created by rja on 07/06/2021.
//

#include <M7_lib/foreach/Foreach.h>
#include "gtest/gtest.h"
#include "M7_lib/basis/AbelianGroup.h"

TEST(AbelianGroup, D2h){
    std::vector<std::string> labels = {"A1g", "B1g", "B2g", "B3g", "A1u", "B1u", "B2u", "B3u"};
    AbelianGroup d2h(labels, [](const size_t& iirrep, const size_t& jirrep){return iirrep^jirrep;});
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
    for (size_t iirrep = A1g; iirrep<=B3u; iirrep++){
        ASSERT_EQ(d2h.product(iirrep, iirrep), A1g);
        ASSERT_EQ(d2h.product(iirrep, A1g), iirrep);
    }
    size_t iirrep = B1g;
    ASSERT_EQ(d2h.product(iirrep, B2g), B3g);
    ASSERT_EQ(d2h.product(iirrep, B3g), B2g);
    ASSERT_EQ(d2h.product(iirrep, A1u), B1u);
    ASSERT_EQ(d2h.product(iirrep, B1u), A1u);
    ASSERT_EQ(d2h.product(iirrep, B2u), B3u);
    ASSERT_EQ(d2h.product(iirrep, B3u), B2u);
}
TEST(AbelianGroup, Loop) {
    std::vector<std::string> labels = {"A1g", "B1g", "B2g", "B3g", "A1u", "B1u", "B2u", "B3u"};
    AbelianGroup d2h(labels, [](const size_t& iirrep, const size_t& jirrep){return iirrep^jirrep;});
    foreach::ctnd::Ordered<2, false, true> cre_pg(d2h.nirrep());
    foreach::ctnd::Ordered<2, false, true> ann_pg(d2h.nirrep());
    foreach::ctnd::Unrestricted<2> cre_spin(2);
    foreach::ctnd::Unrestricted<2> ann_spin(2);
    auto fn = [&](){
        if (d2h.is_conservative(cre_pg.inds(), ann_pg.inds()) && (cre_spin.sum() == ann_spin.sum())) {
            //std::cout << labels[cre_pg[0]] << " " << labels[cre_pg[1]] << " " << labels[ann_pg[0]] << " " << labels[ann_pg[1]] << " " << cre_spin.inds() << " " << ann_spin.inds() << std::endl;
        }
    };
    foreach::ctnd::chain(fn, cre_pg, ann_pg, cre_spin, ann_spin);
}
