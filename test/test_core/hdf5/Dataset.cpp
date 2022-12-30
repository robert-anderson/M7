//
// Created by rja on 07/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/hdf5/IoManager.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/hdf5/File.h"

TEST(Dataset, RealContiguousSaveLoad) {
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
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, mpi::i_am_root());
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read partially on every rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, true);
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read partially on root and last rank if different
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, mpi::i_am_root() || mpi::i_am(mpi::nrank()-1));
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read all on every rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, false, true);
        // no need to gather: all ranks should have the same data
        ASSERT_EQ(save_vec_all, load_vec);
    }
}

TEST(Dataset, ComplexContiguousSaveLoad) {
    const uint_t nitem = hash::in_range(19 + mpi::irank(), 34, 54);
    const uint_t max_nitem_per_op = 7;
    const auto save_vec_real = convert::vector<double>(hash::in_range(123, nitem, 0, 100));
    const auto save_vec_imag = convert::vector<double>(hash::in_range(456, nitem, 0, 100));
    auto save_vec = arith::zip(save_vec_real, save_vec_imag);
    const auto save_vec_all = mpi::all_gatheredv(save_vec);
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_dataset("stuff", save_vec, max_nitem_per_op);
    }
    v_t<std::complex<double>> load_vec;
    {
        /*
         * only reading on the root rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, mpi::i_am_root());
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read partially on every rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, true);
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read partially on root and last rank if different
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, true, mpi::i_am_root() || mpi::i_am(mpi::nrank()-1));
        const auto load_vec_all = mpi::all_gatheredv(save_vec);
        ASSERT_EQ(save_vec_all, load_vec_all);
    }
    load_vec.clear();
    {
        /*
         * read all on every rank
         */
        hdf5::FileReader fr("tmp.h5");
        fr.load_dataset("stuff", load_vec, max_nitem_per_op, false, true);
        // no need to gather: all ranks should have the same data
        ASSERT_EQ(save_vec_all, load_vec);
    }
}