//
// Created by Robert John Anderson on 2020-03-29.
//

#if 0
#include "src/core/basis/DeterminantField.h"
#include "src/core/basis/PermanentField.h"
#include "gtest/gtest.h"
#include "src/core/table/Table.h"
#include "src/core/table/NumericField.h"

struct TestTable1 : public Table {
    NumericField<size_t> test_longs;
    NumericField<float> test_floats;
    NumericField<char> test_chars;
    NumericField<unsigned short> test_shorts;
    //NumericField<std::complex<double>> test_complexes;

    TestTable1(size_t nsegment, size_t nint, size_t nfloat, size_t nchar, size_t nshort, size_t ncomplex) :
        Table("test table", nsegment),
        test_longs(this, nint),
        test_floats(this, nfloat),
        test_chars(this, nchar),
        test_shorts(this, nshort){}//,
        //test_complexes(this, ncomplex) {}
};

TEST(Table, DataIntegrity) {
    const size_t nsegment = 8;
    const size_t nrow = 6;
    const size_t nelement = 9;

    TestTable1 table(nsegment, nelement, nelement, nelement, nelement, nelement);
    table.expand(nrow);
    size_t i = 0;
    for (size_t isegment = 0ul; isegment < nsegment; ++isegment) {
        for (size_t irow = 0ul; irow < nrow; ++irow) {
            for (size_t ielement = 0ul; ielement < nelement; ++ielement) {
                ASSERT_EQ(table.test_longs(irow, isegment, ielement), 0ul);
                ASSERT_EQ(table.test_floats(irow, isegment, ielement), 0);
                ASSERT_EQ(table.test_chars(irow, isegment, ielement), 0);
                ASSERT_EQ(table.test_shorts(irow, isegment, ielement), 0);
                //ASSERT_EQ(table.test_complexes(irow, isegment, ielement), 0);
                table.test_longs(irow, isegment, ielement) = i;
                table.test_floats(irow, isegment, ielement) = i;
                table.test_chars(irow, isegment, ielement) = i; // will overflow
                table.test_shorts(irow, isegment, ielement) = i;
                //table.test_complexes(irow, isegment, ielement) = i;
                ++i;
            }
        }
    }

    i = 0;
    for (size_t isegment = 0ul; isegment < nsegment; ++isegment) {
        for (size_t irow = 0ul; irow < nrow; ++irow) {
            for (size_t ielement = 0ul; ielement < nelement; ++ielement) {
                ASSERT_EQ(table.test_longs(irow, isegment, ielement), i);
                ASSERT_EQ(table.test_floats(irow, isegment, ielement), i);
                ASSERT_EQ(table.test_chars(irow, isegment, ielement), i); // will overflow
                ASSERT_EQ(table.test_shorts(irow, isegment, ielement), i);
                //ASSERT_EQ(table.test_complexes(irow, isegment, ielement), i);
                ++i;
            }
        }
    }
}

struct TestTable2 : public Table {
    DeterminantField test_dets;

    TestTable2(size_t nsegment, size_t nspatorb, size_t ndet) :
        Table("test table", nsegment),
        test_dets(this, ndet, nspatorb) {}
};

TEST(Table, DataIntegrityDeterminants) {
    const size_t nsegment = 8;
    const size_t nrow = 6;
    const size_t nelement = 9;
    const size_t nspatorb = 12;

    TestTable2 table(nsegment, nspatorb, nelement);
    table.expand(nrow);

    std::vector<DeterminantElement> v;
    table.test_dets(0, 0, 0).set(defs::inds{1, 4, 6, 23});
    v.push_back(table.test_dets(0, 0, 0));
    std::cout << v[0].to_string() <<std::endl;
}


struct TestTable3 : public Table {
    DeterminantField test_dets;
    PermanentField test_perms;

    TestTable3(size_t nsegment, size_t nspatorb, size_t ndet, size_t nperm, size_t nmode, size_t nboson_cutoff) :
            Table("test table", nsegment),
            test_dets(this, ndet, nspatorb),
            test_perms(this, nperm, nmode, nboson_cutoff){}
};

TEST(Table, DataIntegrityPermanents) {
    const size_t nsegment = 8;
    const size_t nrow = 6;
    const size_t nelement = 9;
    const size_t nspatorb = 12;
    const size_t nmode = 6;
    const size_t nboson_cutoff = 12;

    TestTable3 table(nsegment, nspatorb, nelement, nelement, nmode, nboson_cutoff);
    table.expand(nrow);

    auto perm_element = table.test_perms(0, 0, 0);
    perm_element(0) = 123;
    perm_element(2) = 44;
    std::cout << utils::to_string(perm_element.to_vector()) << std::endl;
}

#endif