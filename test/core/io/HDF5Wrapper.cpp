//
// Created by rja on 13/12/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/HDF5Wrapper.h"


TEST(HDF5Wrapper, ListWriter){
    hdf5::File f("test.h5", 1);
    auto grp = f.subgroup("container");
    hdf5::ListWriter<float> lw(grp.m_handle, "LIST", {1});
    float v = 123.456;
    lw.write_item(&v);
}

TEST(HDF5Wrapper, Basic) {
#define H5FILE_NAME     "SDS_row.h5"
#define DATASETNAME 	"IntArray"

        /*
         * HDF5 APIs definitions
         */
        hid_t       file_id, dset_id;         /* file and dataset identifiers */
        hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
        int         *data;                    /* pointer to data buffer to write */
        //hsize_t	count[2];	          /* hyperslab selection parameters */
        //hsize_t	offset[2];
        hid_t	plist_id;                 /* property list identifier */
        herr_t	status;

        const hsize_t nelement = 10;
        hsize_t hypslab_count;
        hsize_t hypslab_offset;

        MPI_Comm comm  = MPI_COMM_WORLD;
        MPI_Info info  = MPI_INFO_NULL;

        /*
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);

        /*
         * Create a new file collectively and release property list identifier.
         */
        file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);


        /*
         * Create the dataspace for the dataset.
         */
        filespace = H5Screate_simple(1, &nelement, NULL);

        /*
         * Create the dataset with default properties and close filespace.
         */
        dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /*
         * Each process defines dataset in memory and writes it to the hyperslab
         * in the file.
         */
        hypslab_count = nelement/mpi::nrank();
        hypslab_offset = mpi::irank() * hypslab_count;
        memspace = H5Screate_simple(1, &hypslab_count, NULL);

        ASSERT(hypslab_offset+hypslab_count<=nelement)

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &hypslab_offset, NULL, &hypslab_count, NULL);

        /*
         * Initialize data buffer
         */
        data = (int *) malloc(sizeof(int)*hypslab_count);
        for (size_t i=0; i < hypslab_count; i++) {
            data[i] = (mpi::irank()+1) * 10 + i;
        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
        ASSERT(!status)
        free(data);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
        H5Pclose(plist_id);
        H5Fclose(file_id);

}

TEST(HDF5Wrapper, Complex){
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


TEST(HDF5Wrapper, Vector){
    hdf5::File f("test.h5", 1);
    auto grp = f.subgroup("container");
    hdf5::VectorWriter<float> vw(grp.m_handle, "my_vector", 10, 1, 10);
    std::vector<float> v = {5, 1, 4, 54, 13524, 3, 2134, 43.12, 123.1243, 456};
    vw.write(v.data(), 1);
}
