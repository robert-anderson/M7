//
// Created by jhalson on 30/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/hamiltonian/BosonCouplings.h"
#include "src/core/field/Elements.h"

TEST(BosonCouplings, Element_b0) {
    size_t nboson_cutoff = 4, nmode = 4;
    defs::ham_t v = 0.5, omega = 0.025;

    BosonCouplings bc(nboson_cutoff, nmode, v, omega);

    elements::FermiBosOnv ket(nmode, nmode);
    elements::FermiBosOnv bra(nmode, nmode);

    ket = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    bra = {{1, 2, 3, 4},
           {0, 0, 0, 0}};
    conn::AsFermiBosOnv fbconn(ket, bra);

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
    size_t nboson_cutoff = 4, nmode = 4;
    defs::ham_t v = 0.5, omega = 0.025;

    BosonCouplings bc(nboson_cutoff, nmode, v, omega);

    elements::FermiBosOnv ket(nmode, nmode);
    elements::FermiBosOnv bra(nmode, nmode);

    ket = {{1, 2, 3, 4},
           {2, 4, 0, 1}};
    bra = {{1, 2, 3, 4},
           {2, 4, 0, 2}};
    conn::AsFermiBosOnv fbconn(ket, bra);
    ASSERT_EQ(fbconn.m_bonvconn.nchanged_mode(), 1);
    ASSERT_EQ(fbconn.m_bonvconn.changed_mode(0), 3);
    ASSERT_EQ(fbconn.m_bonvconn.changes(0), 1);

    auto el = bc.get_element_1(fbconn);
    ASSERT_EQ(v, el);
}

#if 0

TEST(BosonCouplings, Element_f1_b1){
    size_t nocc_cutoff=4, nmode=4;
    defs::ham_t v=0.5, omega=0.025;

    BosonCouplings simpleCoupling(nocc_cutoff, nmode, v, omega);

    Determinant dket(nmode);
    Determinant dbra(nmode);

    dket.set(defs::inds{1,2,3,4});
    dbra.set(defs::inds{1,2,3,5});

    AntisymConnection ac(dbra, dket);

    Permanent pket(nmode, nocc_cutoff);
    Permanent pbra(nmode, nocc_cutoff);

    pket = {2,4,0,1};
    pbra = {2,4,0,2};

    PermanentConnection pc(pbra, pket);

    auto el = simpleCoupling.get_element_1(ac, pc);
    ASSERT_EQ(0, el);
}
#endif