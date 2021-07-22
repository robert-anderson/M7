//
// Created by Robert John Anderson on 2020-03-31.
//

#include "gtest/gtest.h"
#include "src/core/basis/FermionOnvConnection.h"
#include "src/core/io/SparseArrayFileReader.h"
#include "src/core/table/BufferedFields.h"
#include "src/core/hamiltonian/Hamiltonian.h"


TEST(Connection, ParticleNumberConserving) {
    const size_t nsite = 70;
    buffered::FermionOnv ket(nsite);
    buffered::FermionOnv bra(nsite);

    defs::inds ketoccorbs = {1, 4, 6, 8, 11, 19, 120, 138, 139};
    defs::inds braoccorbs = {1, 4, 5, 6, 9, 11, 19, 137, 138};
    ASSERT_EQ(ket.m_size, 3 * defs::nbyte_word);
    ASSERT_EQ(ket.m_dsize, 3);
    ASSERT_EQ(ket.m_nbit_in_last_dword, nsite * 2 - 2 * 64);
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
        ASSERT_TRUE(bra == work_det);
        ASSERT_EQ(connection.phase(), value < 0);
    }
}

TEST(Connection, MultiWordPhase) {

    struct FrmConnNew {
        FermionOnvConnection m_conn;
        const size_t m_ndataword;
        std::vector<bool> m_dataword_phases;

        FrmConnNew(size_t nsite) : m_conn(nsite),
                                   m_ndataword(integer_utils::divceil(nsite * 2, defs::nbit_word)),
                                   m_dataword_phases(m_ndataword) {
            m_dataword_phases[0] = false;
        }

        void update_dataword_phases(const FermionOnvField &onv) {
            for (size_t idataword = 1ul; idataword < m_ndataword; ++idataword) {
                auto prev_dataword = onv.get_dataword(idataword - 1);
                bool phase = bit_utils::nsetbit(prev_dataword) & 1ul;
                m_dataword_phases[idataword] = (m_dataword_phases[idataword - 1] != phase);
            }
        }

        /**
         * The individual phase of a bit position within a multi-word bit representation is the number of set bits
         * before that position.
         * @param onv
         * @param ibit
         * @return
         */
        bool independent_phase(const FermionOnvField &onv, const size_t &ibit) {
            auto idataword = ibit / defs::nbit_word;
            auto ibit_in_word = ibit - idataword * defs::nbit_word;
            return m_dataword_phases[idataword] ^
                   (bit_utils::nsetbit_before(onv.get_dataword(idataword), ibit_in_word) & 1ul);
        }

        /**
         * the overall phase of an excitation with respect to the "in" determinant is the product of the independent
         * phases only if the operators are applied in such a way that they do not interfere with one another
         * e.g. if the occupied set is [0, 1, 4, 6, 7, 9]
         * and the excitation is 9 -> 5
         * the independent phase of 9 is true, since there is an odd number of set bits before position 9
         * the independent phase of 5 is true, since there is an odd number of set bits before position 5
         * so overall, the phase is false, which is correct: the electron at position 9 moves past an even number of
         * others to reach its final position
         *
         * on the other hand, if the excitation is 4 -> 8
         * the independent phase of 4 is false, since there is an even number of set bits before position 4
         * the independent phase of 8 is true, since there is an odd number of set bits before position 8
         * so overall, the phase is true, which is clearly incorrect: the electron at position 4 moves past an even
         * number of others to reach its final position
         *
         * in the second example, the number of set positions before 8 is not correct since we have deleted an electron
         * at 4 first, so we need to compute the phase associated with doing the creation operation first.
         *
         * the algorithm is then to work through the creation and annihilation lists, and each time the smallest of the
         * current indices is a creation operator, negate the overall phase if there is an odd number of annihilations
         * remaining
         *
         * @param onv
         * @return
         */
        bool phase(const FermionOnvField &onv) {
            bool out = false;
            update_dataword_phases(onv);
            auto iann = m_conn.nann()-1;
            auto icre = m_conn.ncre()-1;

            while (iann != ~0ul && icre != ~0ul) {
                const auto& cre = m_conn.cre(icre);
                const auto& ann = m_conn.ann(iann);
                if (cre > ann){
                    // next ind is a creation
                    out ^= independent_phase(onv, cre);
                    // take into account the number of annihilation ops jumped-over to apply this op first
                    out ^= !(iann & 1ul);
                    --icre;
                } else {
                    out ^= independent_phase(onv, ann);
                    --iann;
                }
            }
            while (icre != ~0ul) out ^= independent_phase(onv, m_conn.cre(icre--));
            while (iann != ~0ul) out ^= independent_phase(onv, m_conn.ann(iann--));
            /*
             * we worked backwards through the creation operator string, but it is conventionally ordered the other way
             * around, so multiply by the phase of inverting a string of n distinct SQ op indices
             */
            out ^= (m_conn.ncre()/2) & 1;
            return out;
        }
    };

    // choose a large enough nsite so that multiple 64bit datawords are required
    const size_t nsite = 100;
    buffered::FermionOnv bra(nsite);
    buffered::FermionOnv ket(nsite);
    buffered::FermionOnv check_onv(nsite);
    AntisymFermionOnvConnection connection(nsite);
    FrmConnNew frmconn(nsite);


    bra = {{3,  4,  8,  13, 89},
           {13, 78, 95, 98, 99}};

    frmconn.update_dataword_phases(bra);
    ASSERT_FALSE(frmconn.m_wordwise_cum_phase[0]);
    ASSERT_FALSE(frmconn.m_wordwise_cum_phase[1]); // 4 setbits in dataword 0
    ASSERT_FALSE(frmconn.m_wordwise_cum_phase[2]); // 2 setbits in dataword 1
    ASSERT_TRUE(frmconn.m_wordwise_cum_phase[3]);  // 1 setbits in dataword 2


    ASSERT_FALSE(frmconn.independent_phase(bra, 0)); // 0 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 1)); // 0 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 2)); // 0 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 3)); // 0 before
    ASSERT_TRUE(frmconn.independent_phase(bra, 4));  // 1 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 5)); // 2 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 6)); // 2 before

    ASSERT_TRUE(frmconn.independent_phase(bra, nsite));     // 5 before
    ASSERT_TRUE(frmconn.independent_phase(bra, nsite + 13));  // 5 before
    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 14)); // 6 before
    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 15)); // 6 before

    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 77)); // 6 before
    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 78)); // 6 before
    ASSERT_TRUE(frmconn.independent_phase(bra, nsite + 79));  // 7 before

    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 97)); // 8 before
    ASSERT_FALSE(frmconn.independent_phase(bra, nsite + 98)); // 8 before
    ASSERT_TRUE(frmconn.independent_phase(bra, nsite + 99)); // 9 before
    ASSERT_FALSE(frmconn.independent_phase(bra, 2 * nsite)); // 10 before


    ket = bra;
    connection.connect(bra, ket);
    frmconn.m_conn.connect(bra, ket);
    //ASSERT_EQ(connection.phase(), frmconn.phase(bra));
    ASSERT_FALSE(connection.phase());
    connection.apply(bra, check_onv);
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //                              78 -> 3 (-1)         (-1)
    ket = {{3, 4,  8,  13, 89},
           {3, 13, 95, 98, 99}};
    connection.connect(bra, ket);
    frmconn.m_conn.connect(bra, ket);
    ASSERT_EQ(connection.phase(), frmconn.phase(bra));
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), nsite + 78);
    ASSERT_EQ(connection.cre(0), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());

    ASSERT_EQ(ket, check_onv);

    //                              98 -> 3 (-1)         (-1)
    ket = {{3, 4,  8,  13, 89},
           {3, 13, 78, 95, 99}};
    connection.connect(bra, ket);
    frmconn.m_conn.connect(bra, ket);
    ASSERT_EQ(connection.phase(), frmconn.phase(bra));
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), nsite + 98);
    ASSERT_EQ(connection.cre(0), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //     4 -> 99 (-1)             98 -> 3 (-1)         (+1)
    ket = {{3, 8,  13, 89, 99},
           {3, 13, 78, 95, 99}};
    connection.connect(bra, ket);
    frmconn.m_conn.connect(bra, ket);
    ASSERT_EQ(connection.phase(), frmconn.phase(bra));
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(connection.ann(0), 4);
    ASSERT_EQ(connection.ann(1), nsite + 98);
    ASSERT_EQ(connection.cre(0), 99);
    ASSERT_EQ(connection.cre(1), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //     4 -> 99 (-1)             98 -> 14 (+1)        (-1)
    ket = {{3,  8,  13, 89, 99},
           {13, 14, 78, 95, 99}};
    connection.connect(bra, ket);
    frmconn.m_conn.connect(bra, ket);
    ASSERT_EQ(connection.phase(), frmconn.phase(bra));
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), 4);
    ASSERT_EQ(connection.ann(1), nsite + 98);
    ASSERT_EQ(connection.cre(0), 99);
    ASSERT_EQ(connection.cre(1), nsite + 14);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(ket, check_onv);
}
