//
// Created by jhalson on 30/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/hamiltonian/BosonCouplings.h"
#include "src/core/field/Elements.h"

TEST(BosonCouplings, Element_b0) {
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;

    BosonCouplings bc(nboson_cutoff, nsite, v, omega);

    elements::FermiBosOnv ket(nsite);
    elements::FermiBosOnv bra(v);

    ket = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    bra = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    conn::Antisym<1> fbconn(ket, bra);

    auto el = bc.get_element_0(fbconn);
    ASSERT_EQ(el, 0);

    ket = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    bra = {{1, 2, 3, 4},
           {2, 4, 0, 1}};

    fbconn.connect(bra, ket);

    el = bc.get_element_0(fbconn);
    ASSERT_EQ(el, 7 * omega);
}

TEST(BosonCouplings, Element_f0_b1){
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;

    BosonCouplings bc(nboson_cutoff, nsite, v, omega);

    elements::FermiBosOnv ket(nsite);
    elements::FermiBosOnv bra(nsite);

    ket = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    bra = {{1, 2, 3, 4},
           {2, 4, 0, 2}};
    conn::Antisym<1> fbconn(ket, bra);
    ASSERT_EQ(fbconn.m_bonvconn.nchanged_mode(), 1);
    ASSERT_EQ(fbconn.m_bonvconn.changed_mode(0), 3);
    ASSERT_EQ(fbconn.m_bonvconn.changes(0), 1);

    auto el = bc.get_element_1(fbconn);
    ASSERT_EQ(v*std::sqrt(2.0), el);


    ket = {{1, 2, 3, 4},
           {2, 6, 0, 1}};
    bra = {{1, 2, 3, 4},
           {2, 5, 0, 1}};
    fbconn.connect(ket, bra);
    ASSERT_EQ(fbconn.m_bonvconn.nchanged_mode(), 1);
    ASSERT_EQ(fbconn.m_bonvconn.changed_mode(0), 1);
    ASSERT_EQ(fbconn.m_bonvconn.changes(0), -1);

    el = bc.get_element_1(fbconn);
    ASSERT_EQ(v*std::sqrt(6.0), el);
}

TEST(BosonCouplings, Element_f1_b1){
    size_t nboson_cutoff = 4, nsite = 4;
    defs::ham_t v = 0.5, omega = 0.025;

    BosonCouplings bc(nboson_cutoff, nsite, v, omega);

    elements::FermiBosOnv ket(nsite);
    elements::FermiBosOnv bra(nsite);

    ket = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    bra = {{1, 2, 3, 5},
           {2, 4, 0, 2}};
    conn::Antisym<1> fbconn(ket, bra);
    auto el = bc.get_element_1(fbconn);
    ASSERT_EQ(0, el);
}
