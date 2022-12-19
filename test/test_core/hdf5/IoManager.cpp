//
// Created by rja on 07/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/hdf5/IoManager.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/hdf5/File.h"
#include "M7_lib/util/Hash.h"

TEST(IoManager, ContiguousWriteManager) {
    const uint_t nitem = hash::in_range(19 + mpi::irank(), 34, 54);
    const uint_t max_nitem_per_op = 7;
    const auto save_vec = hash::in_range(123, nitem, 0, 100);
    const auto save_vec_all = mpi::all_gatheredv(save_vec);
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_dataset("stuff", save_vec, max_nitem_per_op);
    }
    v_t<hash::digest_t> load_vec;
    {
        /*
         * only reading on the root rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, mpi::i_am_root());
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read on every rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true);
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }

    //hdf5::dataset::VectorSaveManager<int> vsm();
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