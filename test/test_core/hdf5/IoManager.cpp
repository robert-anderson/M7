//
// Created by rja on 07/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/hdf5/IoManager.h"
#include "M7_lib/util/Hash.h"

TEST(IoManager, VectorWriteManager) {
//    const uint_t nitem = 123;
//    const uint_t nitem_per_transfer = 45;
//    const auto ntransfer = integer::divceil(nitem, nitem_per_transfer);
//    const auto nitem_last_transfer = nitem - nitem_per_transfer*(ntransfer-1);
//    const auto src = hash::in_range(0, nitem, 3, 17);
//    hdf5::VectorWriteManager<hash::digest_t> manager(&src, src.size());
//    // allocate buffer to store result of transfers
//    v_t<buf_t> result(nitem * sizeof(hash::digest_t), 0);
//    for (auto d = manager.transfer(); d; d = manager.transfer()) {
//
//    }
}