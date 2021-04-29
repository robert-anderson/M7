//
// Created by rja on 28/04/2021.
//

#include <src/core/hamiltonian/HamiltonianData.h>
#include "gtest/gtest.h"

TEST(HamiltonianData, EncodeDecodeRankLabels){
    size_t ncre = 3, nann = 2;
    auto label = ham_data::encode_rank_label(ncre, nann);
    size_t ncre_chk, nann_chk;
    std::cout << label << std::endl;
    ham_data::decode_rank_label(label, ncre_chk, nann_chk);
    ASSERT_EQ(ncre, ncre_chk);
    ASSERT_EQ(nann, nann_chk);
}