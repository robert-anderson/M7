//
// Created by rja on 13/12/2020.
//

#include <M7_lib/hash/Hashing.h>
#include <M7_lib/field/Row.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/table/Table.h>
#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/table/BufferedFields.h>
#include "gtest/gtest.h"
#include "M7_lib/io/HDF5Wrapper.h"


namespace hdf5_wrapper_test {

}

TEST(HDF5Wrapper, StringVector) {
    auto definitive_irank = hashing::in_range(99, 0, mpi::nrank());
    std::vector<std::string> strings = {"Lorem", "ipsum dolor sit", "amet, consectetur adipiscing", "elit"};
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        gw.save("a_string_vector", strings, definitive_irank);
    }
}

TEST(HDF5Wrapper, String) {
    auto definitive_irank = hashing::in_range(99, 0, mpi::nrank());
    std::string s = "Lorem ipsum dolor sit amet, consectetur adipiscing elit";
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        gw.save("a_string", s, definitive_irank);
    }
}

TEST(HDF5Wrapper, FloatArray) {
    auto definitive_irank = hashing::in_range(99, 0, mpi::nrank());
    defs::inds shape = {2, 3};
    std::vector<float> v = {0.1, 4.5, 1.2, 3, 2.3, 4};
    // make each rank's v differ by one element
    v[2] = hashing::in_range(mpi::irank(), 4, 18);
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        gw.save("a_float_array", v.data(), shape, {"dim0", "dim1"}, definitive_irank);
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        hdf5::GroupReader gr("container", fr);
        auto nelement = nd_utils::nelement(shape);
        ASSERT_EQ(nelement, v.size());
        std::vector<float> v_read(nelement);
        ASSERT_TRUE(gr.child_exists("a_float_array"));
        gr.load("a_float_array", v_read.data(), shape);
        auto v_def = v;
        v_def[2] = hashing::in_range(definitive_irank, 4, 18);
        ASSERT_EQ(v_read, v_def);
    }
}

TEST(HDF5Wrapper, ComplexArray) {
    auto definitive_irank = hashing::in_range(99, 0, mpi::nrank());
    defs::inds shape = {2, 3};
    std::vector<std::complex<float>> v = {{0.1, 1}, {4.5, 2}, {1.2, 3}, {3, 4}, {2.3, 5}, {4, 6}};
    // make each rank's v differ by one element
    v[2].imag(hashing::in_range(mpi::irank(), 4, 18));
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        gw.save("a_complex_array", v.data(), shape, {"dim0", "dim1"}, definitive_irank);
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        hdf5::GroupReader gr("container", fr);
        auto nelement = nd_utils::nelement(shape);
        ASSERT_EQ(nelement, v.size());
        std::vector<std::complex<float>> v_read(nelement);
        gr.load("a_complex_array", v_read.data(), shape);
        auto v_def = v;
        v_def[2].imag(hashing::in_range(definitive_irank, 4, 18));
        ASSERT_EQ(v_read, v_def);
    }
}


//TEST(HDF5Wrapper, BufferedFloats) {
//    auto definitive_irank = hashing::in_range(99, 0, mpi::nrank());
//    buffered::Numbers<float, 2> v({2, 3});
//    v = {1, 2, 3, 4, 5, 6};
//    {
//        hdf5::FileWriter fw("table_test.h5");
//        hdf5::GroupWriter gw("container", fw);
//        gw.save("a_complex_array", v.data(), shape, {"dim0", "dim1"}, definitive_irank);
//    }
//    mpi::barrier();
//    {
//        hdf5::FileReader fr("table_test.h5");
//        hdf5::GroupReader gr("container", fr);
//        auto nelement = nd_utils::nelement(shape);
//        ASSERT_EQ(nelement, v.size());
//        std::vector<std::complex<float>> v_read(nelement);
//        gr.load("a_complex_array", v_read.data(), shape);
//        auto v_def = v;
//        v_def[2].imag(hashing::in_range(definitive_irank, 4, 18));
//        ASSERT_EQ(v_read, v_def);
//    }
//}


TEST(HDF5Wrapper, NumberDistributed) {
    BufferedTable<SingleFieldRow<field::Number<defs::hash_t>>> write_table("test int table", {{"integer_field"}});
    auto read_table = write_table;
    const auto nrow = hashing::in_range(mpi::irank() + 1, 10, 20);
    log::debug_("number of local rows {}", nrow);
    write_table.push_back(nrow);
    auto row = write_table.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        row.m_field = hashing::in_range((row.index() + 1) * (mpi::irank() + 1), 4, 123);
        log::debug_("writing value: {}", defs::hash_t(row.m_field));
    }
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw("container", fw);
        write_table.save(gw, "table");
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        hdf5::GroupReader gr("container", fr);

        ASSERT_EQ(gr.child_name(0), "table");
        ASSERT_EQ(gr.m_parent_handle, fr.m_handle);
        ASSERT_TRUE(gr.child_exists("table"));
        ASSERT_FALSE(gr.child_exists("not_the_table"));
        read_table.load(gr, "table");
    }
    for (row.restart(); row.in_range(); row.step()) {
        ASSERT_EQ(row.m_field, hashing::in_range((row.index() + 1) * (mpi::irank() + 1), 4, 123));
    }
}

//TEST(HDF5Wrapper, Table) {
//
//    struct MyRow : Row {
//        fields::Number<int> m_int;
//        fields::Numbers<double, 3> m_double;
//        fields::FermionOnv m_det;
//
//        MyRow() :
//                m_int(this), m_double(this, {{2,   4,      3},
//                                             {"A", "bbob", "cfs"}}),
//                m_det(this, 9) {}
//    };
//
//    BufferedTable<MyRow> table("Test", {{}});
//    table.push_back(3);
//    table.m_row.restart();
//    table.m_row.m_double = 1.23;
//    table.m_row.m_double[2] = 9.09;
//    table.m_row.step();
//    table.m_row.m_double = 3.55;
//
//    hdf5::FileWriter fw("table_test.h5");
//    hdf5::GroupWriter gw("container", fw);
//    table.write(gw, "table");
//}