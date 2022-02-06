//
// Created by anderson on 2/5/22.
//

#include "HeisenbergFrmHam.h"

HeisenbergFrmHam::HeisenbergFrmHam(defs::ham_t j, Lattice lattice) :
        FrmHam(lattice.nsite(), lattice.nsite(), true), m_j(j), m_lattice(std::move(lattice)){}

HeisenbergFrmHam::HeisenbergFrmHam(const fciqmc_config::FermionHamiltonian &opts) :
        HeisenbergFrmHam(opts.m_heisenberg.m_coupling, lattice::make(opts.m_heisenberg)){}

defs::ham_t HeisenbergFrmHam::get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    /*
     * [ik|jl] only have non-zero coeff if i and k spins are opposite and so also are those of j and l
     * also require site(i)==site(k) and site(j)==site(l)
     */
    if (field::FrmOnv::ispin(i, m_nsite)==field::FrmOnv::ispin(k, m_nsite)) return 0.0;
    if (field::FrmOnv::ispin(j, m_nsite)==field::FrmOnv::ispin(l, m_nsite)) return 0.0;
    auto isite = field::FrmOnv::isite(i, m_nsite);
    auto jsite = field::FrmOnv::isite(j, m_nsite);
    if (isite!=field::FrmOnv::isite(k, m_nsite)) return 0.0;
    if (jsite!=field::FrmOnv::isite(l, m_nsite)) return 0.0;
    // fermi phase not included here
    return 0.5*m_j*m_lattice.m_dense(isite, jsite);
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
            int jspin = onv.get({0, jsite}) ? -1: 1;
            sj_tot+= row.second[ineigh]*jspin;
        }
        int ispin = onv.get({0, isite}) ? -1: 1;
        si_sj_tot+=sj_tot*ispin;
    }
    return si_sj_tot*m_j;
}

defs::ham_t HeisenbergFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
    // fermi phase is always negative
    return -HeisenbergFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
}
