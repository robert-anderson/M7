//
// Created by Robert J. Anderson on 24/07/2021.
//

#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/io/CsvFileReader.h>
#include "M7_lib/connection/FrmOnvConnection.h"
#include "gtest/gtest.h"

namespace frm_onv_connection_test {
    static bool phase_direct(const field::FrmOnv& src, const field::FrmOnv& dst){
        FrmOnvConnection connection(src);
        connection.connect(src, dst);
        return connection.phase(src);
    }
    static bool phase_connect(const field::FrmOnv& src, const field::FrmOnv& dst){
        FrmOnvConnection connection(src);
        FrmOps com(src.m_basis.m_nsite);
        return connection.connect(src, dst, com);
    }
    static bool phase_apply(const field::FrmOnv& src, const field::FrmOnv& dst){
        FrmOnvConnection connection(src);
        connection.connect(src, dst);
        FrmOps com(src.m_basis.m_nsite);
        return connection.apply(src, com);
    }
    static size_t ncre(const field::FrmOnv& src, const field::FrmOnv& dst){
        FrmOnvConnection connection(src);
        connection.connect(src, dst);
        return connection.m_cre.size();
    }
    static size_t nann(const field::FrmOnv& src, const field::FrmOnv& dst){
        FrmOnvConnection connection(src);
        connection.connect(src, dst);
        return connection.m_ann.size();
    }
    static size_t string_chk(const field::FrmOnv& src, const field::FrmOnv& dst, defs::inds_t ann, defs::inds_t cre){
        FrmOnvConnection connection(src);
        connection.connect(src, dst);
        return (connection.m_ann == ann) && (connection.m_cre == cre);
    }
    static size_t string_chk(const field::FrmOnv& src, const field::FrmOnv& dst, defs::inds_t ann, defs::inds_t cre, defs::inds_t com){
        FrmOnvConnection connection(src);
        FrmOps com_chk(src.m_basis.m_nsite);
        connection.connect(src, dst, com_chk);
        return (connection.m_ann == ann) && (connection.m_cre == cre) && (com_chk==com);
    }
    static bool pred_true(bool direct, bool connect, bool apply){
        return direct && connect && apply;
    }
    static bool pred_false(bool direct, bool connect, bool apply){
        return !pred_true(direct, connect, apply);
    }
}

TEST(FrmOnvConnection, Ex01){
    using namespace frm_onv_connection_test;
    const size_t nsite = 5;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {1, 2, 3, 5, 7, 8};
    buffered::FrmOnv dst_onv(nsite);

    dst_onv = {1, 2, 3, 5, 8};
    ASSERT_EQ(ncre(src_onv, dst_onv), 0ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 1ul);
    // even number of electronic exchanges
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));

    dst_onv = {1, 2, 3, 5, 7};
    ASSERT_EQ(ncre(src_onv, dst_onv), 0ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 1ul);
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}


TEST(FrmOnvConnection, Ex10){
    using namespace frm_onv_connection_test;
    const size_t nsite = 5;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {1, 2, 3, 5, 7, 8};
    buffered::FrmOnv dst_onv(nsite);

    dst_onv = {1, 2, 3, 4, 5, 7, 8};
    ASSERT_EQ(ncre(src_onv, dst_onv), 1ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 0ul);
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));

    dst_onv = {1, 2, 3, 5, 7, 8, 9};
    ASSERT_EQ(ncre(src_onv, dst_onv), 1ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 0ul);
    // even number of electronic exchanges
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}


TEST(FrmOnvConnection, Ex11SingleWord){
    using namespace frm_onv_connection_test;
    const size_t nsite = 32;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {{3, 8, 12, 30}, {12, 14, 15, 22}};
    buffered::FrmOnv dst_onv(nsite);
    //                                    12 -> 23 (+)
    dst_onv = {{3, 8, 12, 30}, {14, 15, 22, 23}};
    ASSERT_EQ(ncre(src_onv, dst_onv), 1ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 1ul);
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
    //                                    22 -> 13 (+)
    dst_onv = {{3, 8, 12, 30}, {12, 13, 14, 15}};
    // even number of electronic exchanges
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
    //                                    22 -> 11 (-)
    dst_onv = {{3, 8, 12, 30}, {11, 12, 14, 15}};
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}

TEST(FrmOnvConnection, Ex11MultiWord){
    using namespace frm_onv_connection_test;
    const size_t nsite = 125;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {{3, 8, 12, 30}, {112, 114, 115, 122}};
    buffered::FrmOnv dst_onv(nsite);
    //                                    112 -> 123 (+)
    dst_onv = {{3, 8, 12, 30}, {114, 115, 122, 123}};
    ASSERT_EQ(ncre(src_onv, dst_onv), 1ul);
    ASSERT_EQ(nann(src_onv, dst_onv), 1ul);
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
    //                                   122 -> 113 (+)
    dst_onv = {{3, 8, 12, 30}, {112, 113, 114, 115}};
    // even number of electronic exchanges
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
    //                                    122 -> 111 (-)
    dst_onv = {{3, 8, 12, 30}, {111, 112, 114, 115}};
    // odd number of electronic exchanges
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}

TEST(FrmOnvConnection, Ex22){
    using namespace frm_onv_connection_test;
    const size_t nsite = 123;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {3, 8, 12, 30, 123, 144, 188, 234};
    buffered::FrmOnv dst_onv(nsite);
    dst_onv = {3, 8, 12, 31, 144, 188, 189, 234};
    ASSERT_EQ(src_onv.nsetbit(), 8);
    ASSERT_EQ(dst_onv.nsetbit(), 8);
    ASSERT_EQ(nann(src_onv, dst_onv), 2ul);
    ASSERT_EQ(ncre(src_onv, dst_onv), 2ul);

    ASSERT_TRUE(string_chk(src_onv, dst_onv, {30, 123}, {31, 189}, {3, 8, 12, 144, 188, 234}));

    // an even number of electron coordinates exchanges have occurred
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));

    // an odd number of electron coordinates exchanges have occurred
    dst_onv = {8, 12, 30, 144, 145, 188, 234, 235};
    ASSERT_TRUE(string_chk(src_onv, dst_onv, {3, 123}, {145, 235}, {8, 12, 30, 144, 188, 234}));
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}

TEST(FrmOnvConnection, Ex33){
    using namespace frm_onv_connection_test;
    const size_t nsite = 123;
    buffered::FrmOnv src_onv(nsite);
    src_onv = {3, 8, 12, 30, 123, 144, 188, 234};
    buffered::FrmOnv dst_onv(nsite);
    dst_onv = {3, 8, 12, 32, 123, 144, 190, 235};
    ASSERT_EQ(src_onv.nsetbit(), 8);
    ASSERT_EQ(dst_onv.nsetbit(), 8);
    ASSERT_EQ(nann(src_onv, dst_onv), 3ul);
    ASSERT_EQ(ncre(src_onv, dst_onv), 3ul);

    ASSERT_TRUE(string_chk(src_onv, dst_onv, {30, 188, 234}, {32, 190, 235}, {3, 8, 12, 123, 144}));

    // no electron coordinates have been exchanged
    ASSERT_PRED3(pred_false, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));

    // this time, the parity is true
    dst_onv = {8, 12, 14, 30, 123, 188, 189, 235};
    ASSERT_TRUE(string_chk(src_onv, dst_onv, {3, 144, 234}, {14, 189, 235}, {8, 12, 30, 123, 188}));
    ASSERT_PRED3(pred_true, phase_direct(src_onv, dst_onv), phase_connect(src_onv, dst_onv), phase_apply(src_onv, dst_onv));
}

TEST(Connection, EntireCiPhases) {
    using namespace frm_onv_connection_test;
    /**
     * this file enumerates all determinantal connections in a small CI space along with their associated phases
     */
    NumericCsvFileReader file_reader(defs::c_assets_root + "/parity_test/parity_8.txt", 17);
    defs::inds_t inds(16);
    int value;

    buffered::FrmOnv bra(4);
    buffered::FrmOnv ket(4);
    buffered::FrmOnv work_det(4);

    std::vector<std::string> tokens;
    while (file_reader.next(tokens)) {
        NumericCsvFileReader::parse(tokens.cbegin(), tokens.cbegin()+1, value);
        NumericCsvFileReader::parse(tokens.cbegin()+1, tokens.cend(), inds);
        bra.zero();
        ket.zero();
        for (size_t i = 0ul; i < 8ul; ++i) {
            if (inds[i]) bra.set(i);
        }
        for (size_t i = 8ul; i < 16ul; ++i) {
            if (inds[i]) ket.set(i - 8);
        }
        if (bra.is_zero() || ket.is_zero()) continue;
        if (bra.nsetbit() != ket.nsetbit()) continue;

        if (value<0)
            ASSERT_PRED3(pred_true, phase_direct(ket, bra), phase_connect(ket, bra), phase_apply(ket, bra));
        else
            ASSERT_PRED3(pred_false, phase_direct(ket, bra), phase_connect(ket, bra), phase_apply(ket, bra));
    }
}
