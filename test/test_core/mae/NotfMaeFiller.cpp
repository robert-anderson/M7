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
        { 0,  1,  3,   6, 10, 11}
    };

    /*
     * input the set bit data into the hist table as determinant bitstrings
     */
    for (auto& setbits: setbits_vec) {
        hist.m_row.push_back_jump();
        hist.m_row.m_mbf = setbits;
    }

#if 0
    std::cout << hist.to_string() << std::endl;

    NotfMaeFiller filler(hist);

    auto test_fn = [&](const uintv_t& /*ann_ispinorbs*/, const v_t<uintp_t>& /*ann_isect*/,
                       const uintv_t& /*cre_ispinorbs*/, const v_t<uintp_t>& /*cre_isect*/) -> void {

    };
    filler.fill_foreach_isect_pair(test_fn);
#endif
}