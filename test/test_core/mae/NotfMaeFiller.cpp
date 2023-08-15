//
// Created by anderson on 09/08/2023.
//

#include <test_core/defs.h>
#include <M7_lib/mae/NotfMaeFiller.h>

TEST(NotfMaeFiller, AllIsects) {
    const sys::Basis basis = {{6ul}, {0ul}};
    const NdFormat<c_ndim_wf> wf_fmt({1ul, 1ul});
    buffered::Table<MbfWeightRow> hist("test hist", MbfWeightRow(basis, wf_fmt));

    const v_t<uintv_t> setbits_vec = {
        { 0,  1,  4,   6,  8, 11},
        { 1,  4,  5,   7,  8,  9},
        { 0,  1,  3,   6, 10, 11},
        { 0,  2,  3,   6, 10, 11},
        { 0,  2,  4,   8,  9, 11},
        { 0,  3,  5,   6,  9, 11},
        { 0,  1,  5,   7,  8, 10},
        { 1,  2,  4,   7,  9, 10},
        { 3,  4,  5,   6,  8, 11},
        { 1,  3,  5,   9,  8, 11}
    };

    for (auto& setbits: setbits_vec) {
        hist.m_row.push_back_jump();
        hist.m_row.m_mbf = setbits;
    }
    NotfMaeFiller filler(hist);
    auto test_fn = [&](const uintv_t& ann_ispinorbs, const v_t<uintp_t>& /*ann_isect*/,
                       const uintv_t& cre_ispinorbs, const v_t<uintp_t>& /*cre_isect*/) -> void {
        if (ann_ispinorbs.size()==2) {
            std::cout << ann_ispinorbs[0] << " " << ann_ispinorbs[1] << "   " << cre_ispinorbs[0] << " " << cre_ispinorbs[1] << std::endl;
        }
        else {
            std::cout << ann_ispinorbs[0] << "   " << cre_ispinorbs[0] << std::endl;
        }
    };
    v_t<OpSig> rdm_exsigs;
    rdm_exsigs.emplace_back(opsig::c_doub);
    filler.fill_foreach_set_pair(test_fn, rdm_exsigs);
}


TEST(NotfMaeFiller, HalfExcitPhase) {
    const sys::Basis basis = {{6ul}, {0ul}};
    buffered::Mbf mbf(basis);
    const uintv_t setbits_vec = { 0,  1,  4,   6,  8, 11};

    mbf = setbits_vec;
    conn::Mbf conn(mbf);

    // 0 1 4 6 8 11
    conn.m_ann.add(0);
    ASSERT_EQ(conn.phase(mbf), false);
    // 1 4 6 8 11
    conn.m_ann.add(4);
    ASSERT_EQ(conn.phase(mbf), true);
    // 1 6 8 11
    conn.m_ann.add(6);
    ASSERT_EQ(conn.phase(mbf), false);
    // 1 8 11
    conn.m_ann.add(11);
    ASSERT_EQ(conn.phase(mbf), false);
}
