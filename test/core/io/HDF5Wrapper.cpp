//
// Created by rja on 13/12/2020.
//

#include <src/core/hash/Hashing.h>
#include <src/core/field/Row.h>
#include <src/core/field/Fields.h>
#include <src/core/table/Table.h>
#include <src/core/table/BufferedTable.h>
#include <src/core/table/BufferedFields.h>
#include "gtest/gtest.h"
#include "src/core/io/HDF5Wrapper.h"


namespace hdf5_wrapper_test {

}

TEST(HDF5Wrapper, Number) {
    BufferedTable<SingletRow<fields::Number<int>>> write_table("test int table", {{"integer_field"}});
    auto read_table = write_table;
    const auto nrow = hashing::in_range(mpi::irank()+1, 10, 20);
    write_table.push_back(nrow);
    auto row = write_table.m_row;
    for (row.restart(); row.in_range(); row.step()){
        row.m_field = hashing::in_range((row.m_i+1)*(mpi::irank()+1), 4, 123);
        log::debug_("writing value: {}", (int)row.m_field);
    }
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        write_table.write(gw, "table");
    }

    {
        hdf5::FileReader fr("table_test.h5");
        hdf5::GroupReader gr("container", fr);
        read_table.read(gr, "table");
    }
    for (row.restart(); row.in_range(); row.step()){
        ASSERT_EQ(row.m_field, hashing::in_range((row.m_i+1)*(mpi::irank()+1), 4, 123));
    }

}

TEST(HDF5Wrapper, Table) {

    struct MyRow : Row {
        fields::Number<int> m_int;
        fields::Numbers<double, 3> m_double;
        fields::FermionOnv m_det;

        MyRow():
        m_int(this), m_double(this, {{2, 4, 3}, {"A", "bbob", "cfs"}}),
        m_det(this, 9){}
    };

    BufferedTable<MyRow> table("Test", {{}});
    table.push_back(3);
    table.m_row.restart();
    table.m_row.m_double = 1.23;
    table.m_row.m_double[2] = 9.09;
    table.m_row.step();
    table.m_row.m_double = 3.55;

    hdf5::FileWriter fw("table_test.h5");
    hdf5::GroupWriter gw("container", fw);
    table.write(gw, "table");
}


TEST(HDF5Wrapper, Basic) {
#define H5FILE_NAME     "SDS_row.h5"
#define DATASETNAME    "IntArray"


#if 0

    /*
     * HDF5 APIs definitions
     */
    //hid_t       dset_id;         /* file and dataset identifiers */
    //hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    //hsize_t	count[2];	          /* hyperslab selection parameters */
    //hsize_t	offset[2];
    //herr_t	status;


//    auto get_nelement = [](const std::vector<hsize_t> &v) {
//        hsize_t out = 1;
//        for (auto &i: v) out *= i;
//        return out;
//    };


    typedef int T;
    std::vector<T> data = {9, 9, 9, 9};
    hdf5::FileReader fr(H5FILE_NAME);
    /*

    const size_t nrow = 4;
    defs::inds item_dims_ = {2, 2};

    hdf5::NdListWriter<T> lw(fw, "nd_list_test", item_dims_, nrow);

    std::vector<T> data = {1, 3, 4, 5};

    lw.write_row(2, data.data());
    lw.write_row(0, data.data());
     */
    hdf5::NdListReader<T> lr(fr, "nd_list_test");

    utils::print(data);
    lr.read_row(1, data.data());
    utils::print(data);

    //lw.write(data.data());



    const hsize_t nrow_local = 4;
    const hsize_t nrow_global = nrow_local * mpi::nrank();

    const std::vector<hsize_t> item_dims = {};
    auto list_dims_local = item_dims;
    list_dims_local.insert(list_dims_local.begin(), nrow_local);
    auto list_dims_global = item_dims;
    list_dims_global.insert(list_dims_global.begin(), nrow_global);

    const hsize_t ndim_item = item_dims.size();
    const hsize_t ndim_list = ndim_item + 1;
    //const hsize_t nelement_list_local = get_nelement(list_dims_local);
    //const hsize_t nelement_list_global = get_nelement(list_dims_global);

    std::vector<hsize_t> hypslab_counts(ndim_list, 0ul);
    hypslab_counts = list_dims_local;
    std::vector<hsize_t> hypslab_offsets(ndim_list, 0ul);



    /*
     * Create a new file collectively and release property list identifier.
     */

    hdf5::FileReader fr(H5FILE_NAME);

    hdf5::NdListReader<int> lr (fr, "nd_list_test");
    utils::print(lr.m_list_dims_local);
    utils::print(lr.m_item_dims);

    std::vector<int> data;
    for (size_t i=0; i < nd_utils::nelement(lr.m_list_dims_local); i++) {
        data.push_back((mpi::irank()+1) * 100 + i);
    }

    std::vector<T> read_data(data.size(), 0);
    ASSERT_NE(data, read_data);

    lr.read(read_data.data());
    ASSERT_EQ(data, read_data);


    hdf5::FileWriter fw(H5FILE_NAME);

    const size_t nrow = 4;
    defs::inds item_dims_ = {2, 2};

    hdf5::NdListWriter<T> lw(fw, "nd_list_test", item_dims_, nrow);

    lw.write(data.data());



    ListX list(file.m_handle, "tjtjtj", writemode, item_dims_, nrow);



    if (writemode)
        list.write(data.data());
    else {
        std::vector<T> read_data(data.size(), 0ul);
        list.read(data.data());
        ASSERT_EQ(data, read_data);
    }

    /*
     * Create the dataspace for the dataset.
     */
    filespace = H5Screate_simple(ndim_list, list_dims_global.data(), NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dcreate(group.m_handle, DATASETNAME, H5T_NATIVE_INT, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    hypslab_offsets[0] = mpi::irank() * nrow_local;
    memspace = H5Screate_simple(ndim_list, hypslab_counts.data(), NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, hypslab_offsets.data(), NULL, hypslab_counts.data(), NULL);

    /*
     * Initialize data buffer
     */
    auto dataptr = (int *) malloc(sizeof(int) * nelement_list_local);
    for (size_t i=0; i < nelement_list_local; i++) {
        dataptr[i] = (mpi::irank()+1) * 10 + i;
    }

    /*
     * Create property list for collective dataset write.
     */
    auto plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, dataptr);
    ASSERT(!status)
    free(dataptr);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
#endif
}

TEST(HDF5Wrapper, Complex) {
    std::complex<double> chk = {1.213, 213.346};
    std::complex<double> tmp;

    {
        hdf5::File f("test.h5", 1);
        auto grp = f.subgroup("first level");
        auto grp2 = grp.subgroup("second level");
        grp2.save(chk, "a complex double");
    }
    {
        hdf5::File f("test.h5", 0);
        auto grp = f.subgroup("first level");
        auto grp2 = grp.subgroup("second level");
        grp2.load(tmp, "a complex double");
        ASSERT_FLOAT_EQ(tmp.real(), chk.real());
        ASSERT_FLOAT_EQ(tmp.imag(), chk.imag());
    }
}
