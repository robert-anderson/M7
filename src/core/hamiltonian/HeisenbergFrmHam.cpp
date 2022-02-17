//
// Created by anderson on 2/5/22.
//

#include "HeisenbergFrmHam.h"

HeisenbergFrmHam::HeisenbergFrmHam(defs::ham_t j, Lattice lattice) :
        SpinFrmHam(lattice.nsite(), lattice.nsite(), 0), m_j(j), m_lattice(std::move(lattice)){
    m_contribs_2200.set_nonzero(exsig_utils::ex_double);
    m_contribs_2200.set_nonzero(0);
    log::info("Heisenberg Hamiltonian initialized with J={}; {}", m_j, m_lattice.info());
}

HeisenbergFrmHam::HeisenbergFrmHam(const fciqmc_config::FermionHamiltonian &opts) :
        HeisenbergFrmHam(opts.m_heisenberg.m_coupling, lattice::make(opts.m_heisenberg)){}

defs::ham_t HeisenbergFrmHam::get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    /*
     * normal-ordered SQ operators are:
     *      pa+  qb+  qa   pb
     *  and h.c.:
     *      pb+  qa+  qb   pa
     *
     */
    auto isite = field::FrmOnv::isite(i, m_nsite);
    auto jsite = field::FrmOnv::isite(j, m_nsite);
    // fermi phase not included here, minus sign is due to product of opposite spins
    return -m_j*m_lattice.m_dense(isite, jsite)/2.0;
}

defs::ham_t HeisenbergFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    /*
     * loop over all sites i and accumulate sum of neighboring spins multiplied by the lattice phase factors (in
     * case of periodic BCs - all handled by Lattice class)
     * for each such sum, accumulate its product with the spin of site i into the total energy of the determinant
     * finally, return the accumulation scaled by J
     */
    int si_sj_tot = 0;
    for (size_t isite=0ul; isite<m_nsite; ++isite){
        DEBUG_ASSERT_EQ(onv.site_nocc(isite), 1ul,
                        "spin system is assumed, must not have unoccupied or doubly occupied sites");
        auto row = m_lattice.m_sparse[isite];
        auto nneigh = row.first.size();
        int sj_tot = 0;
        for (size_t ineigh=0ul; ineigh<nneigh; ++ineigh) {
            auto jsite = row.first[ineigh];
            if (jsite<isite) continue; // don't double count contributions
            int jspin = onv.get({0, jsite}) ? 1: -1;
            sj_tot+= row.second[ineigh]*jspin;
        }
        int ispin = onv.get({0, isite}) ? 1: -1;
        si_sj_tot+=sj_tot*ispin;
    }
    // each spin in the product carries a factor of 1/2
    return si_sj_tot*m_j/4.0;
}

defs::ham_t HeisenbergFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
    // fermi phase is always negative
    return -HeisenbergFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[1], conn.m_ann[0]);
}
