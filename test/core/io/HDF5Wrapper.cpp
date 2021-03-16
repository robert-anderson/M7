//
// Created by rja on 13/12/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/HDF5Wrapper.h"


TEST(HDF5Wrapper, ListWriter) {
    hdf5::File f("test.h5", 1);
    auto grp = f.subgroup("container");
    hdf5::ListWriter<float> lw(grp.m_handle, "LIST", {});
    float v = 123.456;
    lw.write_item(&v);
}

TEST(HDF5Wrapper, Basic) {
#define H5FILE_NAME     "SDS_row.h5"
#define DATASETNAME    "IntArray"

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


    typedef int T;
    struct ListX {
        const std::vector<hsize_t> m_item_dims;
        const hsize_t m_ndim_item;
        const hsize_t m_ndim_list;
        const hsize_t m_nrow_local;
        const hsize_t m_nrow_global;
        const std::vector<hsize_t> m_list_dims_local;
        const std::vector<hsize_t> m_list_dims_global;
        const hsize_t m_row_offset;

        std::vector<hsize_t> m_hyperslab_counts;
        std::vector<hsize_t> m_hyperslab_offsets;
        hid_t m_filespace_handle;
        hid_t m_dataset_handle;
        hid_t m_memspace_handle;

        std::vector<hsize_t> get_item_dims(const defs::inds &item_dims) {
            std::vector<hsize_t> out;
            out.reserve(item_dims.size());
            for (auto &i: item_dims) out.push_back(i);
            return out;
        }

        std::vector<hsize_t> get_list_dims_local() {
            std::vector<hsize_t> out;
            out.reserve(m_ndim_list);
            out.push_back(m_nrow_local);
            out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
            return out;
        }

        std::vector<hsize_t> get_list_dims_global() {
            std::vector<hsize_t> out;
            out.reserve(m_ndim_list);
            out.push_back(m_nrow_global);
            out.insert(++out.begin(), m_item_dims.cbegin(), m_item_dims.cend());
            return out;
        }

        hsize_t get_row_offset() {
            std::vector<hsize_t> tmp(mpi::nrank());
            mpi::all_gather(m_nrow_local, tmp);
            hsize_t out = 0ul;
            for (size_t irank = 0ul; irank < mpi::irank(); ++irank) out += tmp[irank];
            return out;
        }

        ListX(hid_t parent, std::string name, bool writemode,
              const defs::inds &item_dims, const size_t &nrow) :
                m_item_dims(get_item_dims(item_dims)),
                m_ndim_item(item_dims.size()),
                m_ndim_list(item_dims.size() + 1),
                m_nrow_local(nrow),
                m_nrow_global(mpi::all_sum(m_nrow_local)),
                m_list_dims_local(get_list_dims_local()),
                m_list_dims_global(get_list_dims_global()),
                m_row_offset(get_row_offset()),
                m_hyperslab_counts(m_ndim_list, 0ul),
                m_hyperslab_offsets(m_ndim_list, 0ul) {
            m_filespace_handle = H5Screate_simple(m_ndim_list, m_list_dims_global.data(), NULL);

            /*
             * Create the dataset with default properties and close filespace.
             */
            m_dataset_handle = H5Dcreate(parent, name.c_str(), hdf5::type<T>(), m_filespace_handle,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(m_filespace_handle);

            /*
             * Each process defines dataset in memory and writes it to the hyperslab
             * in the file.
             */
            m_hyperslab_offsets[0] = m_row_offset;
            m_hyperslab_counts = m_list_dims_local;
            m_memspace_handle = H5Screate_simple(m_ndim_list, m_hyperslab_counts.data(), NULL);


            m_filespace_handle = H5Dget_space(m_dataset_handle);
            H5Sselect_hyperslab(m_filespace_handle, H5S_SELECT_SET, m_hyperslab_offsets.data(),
                                NULL, m_hyperslab_counts.data(), NULL);
        }


        ~ListX() {
            H5Sclose(m_filespace_handle);
            H5Dclose(m_dataset_handle);
            H5Sclose(m_memspace_handle);
        }

        void write(const T *data) {
            auto plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

            auto status = H5Dwrite(m_dataset_handle, hdf5::type<T>(), m_memspace_handle,
                                   m_filespace_handle, plist_id, data);
            ASSERT(!status)
            H5Pclose(plist_id);
        }

        void read(T *data) {
            auto plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

            auto status = H5Dread(m_dataset_handle, hdf5::type<T>(), m_memspace_handle,
                                  m_filespace_handle, plist_id, data);
            ASSERT(!status)
            H5Pclose(plist_id);
        }

    };



    /*
     * Create a new file collectively and release property list identifier.
     */

    hdf5::FileReader fr(H5FILE_NAME);

    hdf5::NdListReader<T> lr (fr, "nd_list_test");
    utils::print(lr.m_list_dims_local);
    utils::print(lr.m_item_dims);

    std::vector<T> data;
    for (size_t i=0; i < nd_utils::nelement(lr.m_list_dims_local); i++) {
        data.push_back((mpi::irank()+1) * 100 + i);
    }

    std::vector<T> read_data(data.size(), 0);
    ASSERT_NE(data, read_data);

    lr.read(read_data.data());
    ASSERT_EQ(data, read_data);


#if 0
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


TEST(HDF5Wrapper, Vector) {
    hdf5::File f("test.h5", 1);
    auto grp = f.subgroup("container");
    hdf5::VectorWriter<float> vw(grp.m_handle, "my_vector", 10, 1, 10);
    std::vector<float> v = {5, 1, 4, 54, 13524, 3, 2134, 43.12, 123.1243, 456};
    vw.write(v.data(), 1);
}
