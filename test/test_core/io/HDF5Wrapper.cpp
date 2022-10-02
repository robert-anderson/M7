//
// Created by Robert J. Anderson on 13/12/2020.
//

#include "gtest/gtest.h"

#include <M7_lib/util/Hash.h>
#include <M7_lib/field/Row.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/hdf5/Dataset.h>

TEST(HDF5Wrapper, NativeTypes) {
    using namespace hdf5;
    /*
     * char
     * signed char
     * unsigned char
     * short
     * unsigned short
     * int
     * unsigned int
     * long
     * unsigned long
     * long long
     * unsigned long long
     * float
     * double
     */
    /*
     * first, the native types
     */
    ASSERT_EQ(type_ind<char>(), 1ul);
    ASSERT_EQ(type_ind<signed char>(), 2ul);
    ASSERT_EQ(type_ind<unsigned char>(), 3ul);
    ASSERT_EQ(type_ind<short>(), 4ul);
    ASSERT_EQ(type_ind<unsigned short>(), 5ul);
    ASSERT_EQ(type_ind<int>(), 6ul);
    ASSERT_EQ(type_ind<unsigned int>(), 7ul);
    ASSERT_EQ(type_ind<long>(), 8ul);
    ASSERT_EQ(type_ind<unsigned long>(), 9ul);
    ASSERT_EQ(type_ind<long long>(), 10ul);
    ASSERT_EQ(type_ind<unsigned long long>(), 11ul);
    ASSERT_EQ(type_ind<float>(), 12ul);
    ASSERT_EQ(type_ind<double>(), 13ul);
    /*
     * now we don't know which native types correspond to the fixed-width integer types on any given machine, so just
     * assert that their type indices are non-zero
     */
    ASSERT_TRUE(type_ind<int8_t>());
    ASSERT_TRUE(type_ind<int16_t>());
    ASSERT_TRUE(type_ind<int32_t>());
    ASSERT_TRUE(type_ind<int64_t>());
    ASSERT_TRUE(type_ind<uint8_t>());
    ASSERT_TRUE(type_ind<uint16_t>());
    ASSERT_TRUE(type_ind<uint32_t>());
    ASSERT_TRUE(type_ind<uint64_t>());
}

TEST(HDF5Wrapper, StringType) {
    str_t s = "Lorem ipsum dolor sit";
    hdf5::Type type(&s);
    ASSERT_EQ(type.m_size, s.size());
}

//TEST(HDF5Wrapper, StringVector) {
//    auto definitive_irank = hash::in_range(99, 0, mpi::nrank());
//    strv_t strings = {"Lorem", "ipsum dolor sit", "amet, consectetur adipiscing", "elit"};
//    {
//        hdf5::FileWriter fw("table_test.h5");
//        hdf5::GroupWriter gw(fw, "container");
//        gw.write_data("a_string_vector", {1}, hdf5::Type(&strings));
//        dw.write(strings);
//    }
//}

TEST(HDF5Wrapper, String) {
    auto definitive_irank = hash::in_range(99, 0, mpi::nrank());
    str_t s = "Lorem ipsum dolor sit amet, consectetur adipiscing elit";
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw(fw, "container");
        gw.write_data("a_string", s, definitive_irank);
    }
}

TEST(HDF5Wrapper, FloatArray) {
    auto definitive_irank = hash::in_range(99, 0, mpi::nrank());
    uintv_t shape = {2, 3};
    v_t<float> v = {0.1, 4.5, 1.2, 3, 2.3, 4};
    // make each rank's v differ by one element
    v[2] = hash::in_range(mpi::irank(), 4, 18);
    {
        hdf5::FileWriter fw("table_test.h5");
        fw.write_data("a_float_array", v);
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        auto nelement = nd::nelement(shape);
        ASSERT_EQ(nelement, v.size());
        ASSERT_TRUE(fr.child_exists("a_float_array"));
        auto v_read = fr.read_data<v_t<float>>("a_float_array");
        auto v_def = v;
        v_def[2] = hash::in_range(definitive_irank, 4, 18);
        ASSERT_EQ(v_read, v_def);
    }
}

TEST(HDF5Wrapper, ComplexArray) {
    auto definitive_irank = hash::in_range(99, 0, mpi::nrank());
    uintv_t shape = {2, 3};
    v_t<std::complex<float>> v = {{0.1, 1}, {4.5, 2}, {1.2, 3}, {3, 4}, {2.3, 5}, {4, 6}};
    // make each rank's v differ by one element
    v[2].imag(hash::in_range(mpi::irank(), 4, 18));
    {
        hdf5::FileWriter fw("table_test.h5");
        fw.write_data("a_complex_array", v.data(), shape, {"dim0", "dim1"}, definitive_irank);
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        const auto v_read = fr.read_data<v_t<std::complex<float>>>("a_complex_array");
        auto v_def = v;
        v_def[2].imag(hash::in_range(definitive_irank, 4, 18));
        ASSERT_EQ(v_read, v_def);
    }
}


//TEST(HDF5Wrapper, BufferedFloats) {
//    auto definitive_irank = hash::in_range(99, 0, mpi::nrank());
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
//        v_t<std::complex<float>> v_read(nelement);
//        gr.load("a_complex_array", v_read.data(), shape);
//        auto v_def = v;
//        v_def[2].imag(hash::in_range(definitive_irank, 4, 18));
//        ASSERT_EQ(v_read, v_def);
//    }
//}


TEST(HDF5Wrapper, NumberDistributed) {
    buffered::Table<SingleFieldRow<field::Number<hash::digest_t>>> write_table("test int table", {"integer_field"});
    ASSERT_EQ(write_table.nrecord(), 0ul);
    auto read_table = write_table;
    const auto nrow = hash::in_range(mpi::irank() + 1, 10, 20);
    logging::debug_("number of local rows {}", nrow);
    write_table.push_back(nrow);
    auto row = write_table.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        row.m_field = hash::in_range((row.index() + 1) * (mpi::irank() + 1), 4, 123);
        logging::debug_("writing value: {}", hash::digest_t(row.m_field));
    }
    {
        hdf5::FileWriter fw("table_test.h5");
        hdf5::GroupWriter gw(fw, "container");
        write_table.save(gw, "table");
    }
    mpi::barrier();
    {
        hdf5::FileReader fr("table_test.h5");
        hdf5::GroupReader gr(fr, "container");

        ASSERT_EQ(gr.child_name(0), "table");
        ASSERT_TRUE(gr.child_exists("table"));
        ASSERT_FALSE(gr.child_exists("not_the_table"));
        read_table.load(gr, "table");
    }
    for (row.restart(); row.in_range(); row.step()) {
        ASSERT_EQ(row.m_field, hash::in_range((row.index() + 1) * (mpi::irank() + 1), 4, 123));
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