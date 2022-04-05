//
// Created by anderson on 2/8/22.
//

#include <M7_lib/field/Mbf.h>
#include "gtest/gtest.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"

namespace heisenberg_test {
    static std::vector<defs::inds> spinvecs() {
        return {{1, 1, 1, 0, 0, 0}, {1, 1, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 1},
                {1, 0, 1, 1, 0, 0}, {1, 0, 1, 0, 1, 0}, {1, 0, 1, 0, 0, 1}, {1, 0, 0, 1, 1, 0},
                {1, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 1, 1, 0, 1, 0},
                {0, 1, 1, 0, 0, 1}, {0, 1, 0, 1, 1, 0}, {0, 1, 0, 1, 0, 1}, {0, 1, 0, 0, 1, 1},
                {0, 0, 1, 1, 1, 0}, {0, 0, 1, 1, 0, 1}, {0, 0, 1, 0, 1, 1}, {0, 0, 0, 1, 1, 1}};
    }

    static std::vector<defs::ham_comp_t> energies() {
        return {0.5, -0.5, -0.5,  0.5, -0.5, -1.5, -0.5, -0.5, -0.5,  0.5,  0.5,
                -0.5, -0.5, -0.5, -1.5, -0.5,  0.5, -0.5, -0.5,  0.5};
    }

    static void set_onv_from_spinvec(field::FrmOnv& onv, const defs::inds& spinvec){
        onv.zero();
        size_t isite = 0ul;
        for (auto &spin: spinvec) onv.set({spin, isite++});
    }
}

TEST(HeisenbergFrmHam, LocalExchangeOnly){
    Lattice::Spec spec(Lattice::Ortho, {6}, {1});
    OrthoLattice lattice(spec);
    for (size_t irow=0ul; irow<lattice.nsite(); ++irow) {
        ASSERT_EQ(lattice.m_sparse.nentry(irow), 2);
    }
    HeisenbergFrmHam ham(1, lattice);
    buffered::FrmOnv src(ham.m_nsite);
    buffered::FrmOnv dst(ham.m_nsite);
    conn::FrmOnv conn(src);

    using namespace heisenberg_test;

    set_onv_from_spinvec(src, {1, 0, 0, 1, 1, 0});
    set_onv_from_spinvec(dst, {1, 0, 0, 1, 0, 1});
    conn.connect(src, dst);
    ASSERT_TRUE(ham.get_element_2200(src, conn));

    set_onv_from_spinvec(dst, {0, 0, 0, 1, 1, 1});
    conn.connect(src, dst);
    // with PBCs, boundary sites on the lattice with opposite spins are allowed to exchange
    ASSERT_TRUE(ham.get_element_2200(src, conn));

    set_onv_from_spinvec(dst, {1, 1, 0, 1, 0, 0});
    conn.connect(src, dst);
    // non-neighboring sites on the lattice with opposite spins are not allowed to exchange
    ASSERT_FALSE(ham.get_element_2200(src, conn));
}

TEST(HeisenbergFrmHam, Elements){
    Lattice::Spec spec(Lattice::Ortho, {6}, {1});
    OrthoLattice lattice(spec);
    for (size_t irow=0ul; irow<lattice.nsite(); ++irow) {
        ASSERT_EQ(lattice.m_sparse.nentry(irow), 2);
    }
    HeisenbergFrmHam ham(1, lattice);
    buffered::FrmOnv src(ham.m_nsite);
    buffered::FrmOnv dst(ham.m_nsite);
    conn::FrmOnv conn(src);

    auto spinvecs = heisenberg_test::spinvecs();
    auto energies = heisenberg_test::energies();

    size_t noffdiag_nonzero = 0ul;
    for (size_t isrc=0ul; isrc<spinvecs.size(); ++isrc) {
        heisenberg_test::set_onv_from_spinvec(src, spinvecs[isrc]);
        ASSERT_FLOAT_EQ(ham.get_energy(src), energies[isrc]);
        for (auto & dst_spinvec : spinvecs) {
            heisenberg_test::set_onv_from_spinvec(dst, dst_spinvec);
            conn.connect(src, dst);
            if (conn.exsig()!=exsig_utils::ex_double) continue;
            auto helem = ham.get_element_2200(src, conn);
            if (!helem) continue;
            ASSERT_FLOAT_EQ(helem, 0.5);
            ++noffdiag_nonzero;
        }
    }
    /*
     * check that we have the expected number of nonzero offdiagonal elements
     */
    ASSERT_EQ(noffdiag_nonzero, 72ul);
}