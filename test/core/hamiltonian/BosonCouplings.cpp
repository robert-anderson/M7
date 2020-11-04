//
// Created by jhalson on 30/09/2020.
//

#if 0

#include "src/core/hamiltonian/BosonCouplings.h"
#include "gtest/gtest.h"
#include "src/core/basis/PermanentConnection.h"
#include "src/core/basis/Permanent.h"

TEST(BosonCouplings, Element_b0){
    size_t nocc_cutoff=4, nmode=4;
    defs::ham_t V=0.5, omega=0.025;

    BosonCouplings simpleCoupling(nocc_cutoff, nmode, V, omega);

    Permanent ket(nmode, nocc_cutoff);
    Permanent bra(nmode, nocc_cutoff);
    PermanentConnection pc(bra, ket);

    auto el = simpleCoupling.get_element_0(pc);
    ASSERT_EQ(el, 0);

    ket = {2,4,0,1};
    bra = {2,4,0,1};

    pc.connect(bra, ket);

    el = simpleCoupling.get_element_0(pc);
    ASSERT_EQ(el, 7*omega);
}

TEST(BosonCouplings, Element_f0_b1){
    size_t nocc_cutoff=4, nmode=4;
    defs::ham_t V=0.5, omega=0.025;

    BosonCouplings simpleCoupling(nocc_cutoff, nmode, V, omega);

    Determinant dket(nmode);

    dket.set(defs::inds{1,2,3,4});
    Determinant dbra = dket;

    AntisymConnection ac(dbra, dket);

    Permanent pket(nmode, nocc_cutoff);
    Permanent pbra(nmode, nocc_cutoff);

    pket = {2,4,0,1};
    pbra = {2,4,0,2};

    PermanentConnection pc(pbra, pket);

    auto el = simpleCoupling.get_element_1(ac, pc);
    ASSERT_EQ(V*ac.ncom(), el);
}

TEST(BosonCouplings, Element_f1_b1){
    size_t nocc_cutoff=4, nmode=4;
    defs::ham_t V=0.5, omega=0.025;

    BosonCouplings simpleCoupling(nocc_cutoff, nmode, V, omega);

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
    ASSERT_EQ(V, el);
}

#endif