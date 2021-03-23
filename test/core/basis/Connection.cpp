//
// Created by Robert John Anderson on 2020-03-31.
//

#include "gtest/gtest.h"
#include "src/core/basis/FermionOnvConnection.h"
#include "src/core/io/SparseArrayFileReader.h"
#include "src/core/table/BufferedFields.h"

TEST(Connection, ParticleNumberConserving){
    const size_t nsite = 70;
    buffered::FermionOnv ket(nsite);
    buffered::FermionOnv bra(nsite);

    defs::inds ketoccorbs = {1, 4, 6, 8, 11, 19, 120, 138, 139};
    defs::inds braoccorbs = {1, 4, 5, 6, 9, 11, 19, 137, 138};
    ASSERT_EQ(ket.m_size, 3*defs::nbyte_data);
    ASSERT_EQ(ket.m_item_dsize, 3);
    ASSERT_EQ(ket.m_nbit_in_last_dword, nsite*2 - 2*64);
    ket = ketoccorbs;
    bra = braoccorbs;

    ASSERT_EQ(ket.nsetbit(), 9);
    ASSERT_EQ(bra.nsetbit(), 9);

    FermionOnvConnection conn(ket, bra);

    ASSERT_EQ(conn.nann(), 3);
    ASSERT_EQ(conn.ncre(), 3);

    ASSERT_EQ(conn.ann(0), 8);
    ASSERT_EQ(conn.ann(1), 120);
    ASSERT_EQ(conn.ann(2), 139);

    ASSERT_EQ(conn.cre(0), 5);
    ASSERT_EQ(conn.cre(1), 9);
    ASSERT_EQ(conn.cre(2), 137);

    AntisymFermionOnvConnection aconn(ket, bra);
    ASSERT_EQ(aconn.ncom(), 6);
    ASSERT_EQ(aconn.com(0), 1);
    ASSERT_EQ(aconn.com(1), 4);
    ASSERT_EQ(aconn.com(2), 6);
    ASSERT_EQ(aconn.com(3), 11);
    ASSERT_EQ(aconn.com(4), 19);
    ASSERT_EQ(aconn.com(5), 138);

    ASSERT_FALSE(aconn.phase());
}

TEST(Connection, Phase) {

    SparseArrayFileReader<float> file_reader(
            defs::assets_root + "/parity_test/parity_8.txt",
            16ul, false, false);
    defs::inds inds(16);
    float value;

    buffered::FermionOnv bra(4);
    buffered::FermionOnv ket(4);
    buffered::FermionOnv work_det(4);
    AntisymFermionOnvConnection connection(ket);

    while (file_reader.next(inds, value)) {
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
        connection.connect(ket, bra);
        ASSERT_EQ(connection.phase(), value < 0);
        connection.apply(ket, work_det);
        ASSERT_TRUE(bra==work_det);
        ASSERT_EQ(connection.phase(), value < 0);
    }
}


TEST(Connection, New){

    using namespace fields;
    struct Onv : MultiField<FermionOnv, BosonOnv> {
        FermionOnv &m_fonv;
        BosonOnv &m_bonv;

        Onv(Row *row, size_t nsite) :
                MultiField<FermionOnv, BosonOnv>(row, {nullptr, nsite}, {nullptr, nsite}),
                m_fonv(get<0>()), m_bonv(get<1>()) {
        }

        Onv(const Onv& other): Onv(other.m_fonv.row_of_copy(), other.m_fonv.m_nsite){}

        Onv &operator=(const Onv& other) {
            m_fonv = other.m_fonv;
            m_bonv = other.m_bonv;
            return *this;
        }

        Onv &operator=(const std::pair<defs::inds, defs::inds>& inds) {
            m_fonv = inds.first;
            m_bonv = inds.second;
            return *this;
        }
    };

    struct Onvs : MultiField<FermionOnvs, BosonOnvs> {
        FermionOnvs &m_fonv;
        BosonOnvs &m_bonv;

        Onvs(Row *row, size_t nitem, size_t nsite) :
                MultiField<FermionOnvs, BosonOnvs>(row, {nullptr, nitem, nsite}, {nullptr, nitem, nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)) {}
    };


    struct FermionConnection {
        defs::inds m_ann;
        defs::inds m_cre;
        defs::inds m_com;
        bool m_phase;

        FermionConnection(size_t nsite){
            m_ann.reserve(2*nsite);
            m_cre.reserve(2*nsite);
            m_com.reserve(2*nsite);
        }

        size_t nann() const {
            return m_ann.size();
        }
        size_t ncre() const {
            return m_cre.size();
        }
        size_t ncom() const {
            return m_com.size();
        }

        void zero(){
            m_cre.clear();
            m_ann.clear();
            m_com.clear();
        }

        void add_cre(const size_t &i){
            ASSERT(m_cre.size()<m_cre.capacity());
            m_cre.push_back(i);
        }

        void add_ann(const size_t &i){
            ASSERT(m_ann.size()<m_ann.capacity());
            m_ann.push_back(i);
        }

        void add(const size_t &ann, const size_t &cre){
            add_ann(ann);
            add_cre(cre);
        }

        void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2){
            add_ann(ann1);
            add_ann(ann2);
            add_cre(cre1);
            add_cre(cre2);
        }

        void sort(){
            std::sort(m_cre.begin(), m_cre.end());
            std::sort(m_ann.begin(), m_ann.end());
        }

        void connect(const FermionOnv &in, const FermionOnv &out) {
            ASSERT(!in.is_zero())
            zero();

            defs::data_t in_work, out_work, work;
            for (size_t idataword = 0ul; idataword<in.m_item_dsize; ++idataword){
                in_work = in.get_dataword(idataword);
                out_work = out.get_dataword(idataword);
                work = in_work &~ out_work;
                while (work) add_ann(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
                work = out_work &~ in_work;
                while (work) add_cre(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
            }


            size_t nperm = 0ul;
            auto ann_iter = m_ann.begin();
            auto cre_iter = m_cre.begin();

            for (size_t idataword = 0ul; idataword<in.m_item_dsize; ++idataword){
                in_work = in.get_dataword(idataword);
                out_work = out.get_dataword(idataword);
                work = in_work & out_work;
                while (work) {
                    auto setbit = bit_utils::next_setbit(work) + idataword * defs::nbit_data;
                    while (ann_iter != m_ann.end() && *ann_iter < setbit) {
                        // an annihilation operator has been passed in the iteration over common indices
                        ann_iter++;
                        nperm += ncom();
                    }
                    while (cre_iter != m_cre.end() && *cre_iter < setbit) {
                        // a creation operator has been passed in the iteration over common indices
                        cre_iter++;
                        nperm += ncom();
                    }
                    m_com.push_back(setbit);
                }
            }
            while (ann_iter != m_ann.end()) {ann_iter++; nperm += ncom();}
            while (cre_iter != m_cre.end()) {cre_iter++; nperm += ncom();}
            m_phase = nperm & 1ul;
        }
    };

    struct BosonConnection {

    };

    struct ConnectionGeneral {
        FermionConnection m_fermion;
        BosonConnection m_boson;
    };

    struct Promoter {

    };

//    struct PromoterSet {
//        ConnectionGeneral m_promoted_conn;
//        std::vector<Promoter> m_proms_f;
//        PromoterSet(size_t max_ncom_f, size_t max_ncom_b){
//
//        }
//    };



    SparseArrayFileReader<float> file_reader(
            defs::assets_root + "/parity_test/parity_8.txt",
            16ul, false, false);
    defs::inds inds(16);
    float value;

    buffered::FermionOnv bra(4);
    buffered::FermionOnv ket(4);
    buffered::FermionOnv work_det(4);
    FermionConnection connection(4);

    while (file_reader.next(inds, value)) {
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
        connection.connect(ket, bra);
        ASSERT_EQ(connection.m_phase, value < 0);
//        connection.apply(ket, work_det);
//        ASSERT_TRUE(bra==work_det);
//        ASSERT_EQ(connection.phase(), value < 0);
    }
}