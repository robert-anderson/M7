//
// Created by anderson on 2/8/22.
//

#include <src/core/field/Mbf.h>
#include "gtest/gtest.h"
#include "src/core/hamiltonian/HeisenbergFrmHam.h"

namespace heisenberg_test {
    static std::vector<defs::inds> spinvecs() {
        return {{1, 1, 1, 0, 0, 0}, {1, 1, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 1},
                {1, 0, 1, 1, 0, 0}, {1, 0, 1, 0, 1, 0}, {1, 0, 1, 0, 0, 1}, {1, 0, 0, 1, 1, 0},
                {1, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 1, 1, 0, 1, 0},
                {0, 1, 1, 0, 0, 1}, {0, 1, 0, 1, 1, 0}, {0, 1, 0, 1, 0, 1}, {0, 1, 0, 0, 1, 1},
                {0, 0, 1, 1, 1, 0}, {0, 0, 1, 1, 0, 1}, {0, 0, 1, 0, 1, 1}, {0, 0, 0, 1, 1, 1}};
    }

    static std::vector<defs::ham_comp_t> energies() {
        return {1.0, -1.0, -1.0, 1.0, -1.0, -3.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -3.0, -1.0, 1.0, -1.0,
                -1.0, 1.0};
    }

    static void set_onv_from_spinvec(field::FrmOnv& onv, const defs::inds& spinvec){
        onv.zero();
        size_t isite = 0ul;
        for (auto &spin: spinvec) onv.set({spin, isite++});
    }
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
            ASSERT_FLOAT_EQ(helem, 1.0);
            ++noffdiag_nonzero;
        }
    }
    /*
     * check that we have the expected number of nonzero offdiagonal elements
     */
    ASSERT_EQ(noffdiag_nonzero, 72ul);
}