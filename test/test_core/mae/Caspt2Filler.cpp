//
// Created by Robert John Anderson on 19/10/2023.
//

#include "test_core/defs.h"
#include "M7_lib/mae/Caspt2Filler.h"

TEST(Caspt2Filler, Rdm1) {
    conf::Document doc;
    doc.m_av_ests.m_rdm.m_ranks = {"1"};
    const uint_t nelec = 6ul;
    const sys::frm::Electrons electrons(nelec);
    const sys::Particles particles {electrons, sys::bos::Bosons(0ul, true)};
    const sys::Basis basis = {{6ul}, {0ul}};
    const sys::Sector sector(basis, particles);

    const NdFormat<c_ndim_wf> wf_fmt({1ul, 1ul});
    buffered::Table<MbfWeightRow> hist("test hist", MbfWeightRow(basis, wf_fmt));

    // todo: add more weights
    const v_t<std::pair<wf_t, uintv_t>> weights_setbits_vec = {
            {1.0, { 0,  1,  4,   6,  8, 11}},
            {1.3, { 1,  4,  5,   7,  8,  9}},
            {-0.2, { 0,  1,  3,   6, 10, 11}},
            { 0,  2,  4,   8,  9, 11},
            { 0,  3,  5,   6,  9, 11},
            { 0,  1,  5,   7,  8, 10},
            { 1,  2,  4,   7,  9, 10},
            { 3,  4,  5,   6,  8, 11},
            { 1,  3,  5,   9,  8, 11}
    };

    for (auto& pair: weights_setbits_vec) {
        hist.m_row.push_back_jump();
        hist.m_row.m_weight = pair.first;
        hist.m_row.m_mbf = pair.second;
    }

    Caspt2Filler filler(hist, nullptr);

    PureRdm rdm1(doc.m_av_ests.m_rdm, opsig::c_sing, sector, 1, "1RDM_test") ;

    filler.fill_rdm1(&rdm1);
}