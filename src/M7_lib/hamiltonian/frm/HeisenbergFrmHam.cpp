//
// Created by Robert J. Anderson on 2/5/22.
//

#include "HeisenbergFrmHam.h"


HeisenbergFrmHam::HeisenbergFrmHam(defs::ham_t j, const std::shared_ptr<lattice::Base> &lattice) :
    SpinModelFrmHam(lattice), m_j(j){
    m_contribs_2200.set_nonzero(utils::exsig::ex_double);
    m_contribs_2200.set_nonzero(0);
    log::info("Heisenberg Hamiltonian initialized with J={}; {}", m_j, m_basis.m_lattice->info());
}


HeisenbergFrmHam::HeisenbergFrmHam(FrmHam::opt_pair_t opts) :
        HeisenbergFrmHam(opts.m_ham.m_heisenberg.m_coupling, lattice::make(opts.m_ham.m_heisenberg)){}

defs::ham_t HeisenbergFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    /*
     * normal-ordered SQ operators are:
     *      pa+  qb+  qa   pb
     *  and h.c.:
     *      pb+  qa+  qb   pa
     *
     */
    const auto asite = m_basis.isite(a);
    const auto bsite = m_basis.isite(b);
    const auto isite = m_basis.isite(i);
    const auto jsite = m_basis.isite(j);
    if (asite == isite){
        if (bsite != jsite) return 0.0;
    }
    else if (bsite == isite) {
        if (asite != jsite) return 0.0;
    }
    else {
        return 0.0;
    }
    // fermi phase not included here, minus sign is due to product of opposite spins
    return -m_j * m_basis.m_lattice->phase(asite, bsite) / 2.0;
}

defs::ham_t HeisenbergFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    /*
     * loop over all sites i and accumulate sum of neighboring spins multiplied by the lattice phase factors (in
     * case of periodic BCs - all handled by Lattice class)
     * for each such sum, accumulate its product with the spin of site i into the total energy of the determinant
     * finally, return the accumulation scaled by J
     */
    int si_sj_tot = 0;
    for (size_t isite=0ul; isite<m_basis.m_nsite; ++isite){
        DEBUG_ASSERT_EQ(onv.site_nocc(isite), 1ul,
                        "spin system is assumed, must not have unoccupied or doubly occupied sites");
        m_basis.m_lattice->get_adj_row(isite, m_work_adj_row);
        int sj_tot = 0;
        for (const auto& elem: m_work_adj_row){
            auto jsite = elem.m_isite;
            if (jsite<isite) continue; // don't double count contributions
            int jspin = onv.get({0, jsite}) ? 1: -1;
            sj_tot+= elem.m_phase*jspin;
        }
        int ispin = onv.get({0, isite}) ? 1: -1;
        si_sj_tot+=sj_tot*ispin;
    }
    // each spin in the product carries a factor of 1/2
    return si_sj_tot*m_j/4.0;
}

defs::ham_t HeisenbergFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), utils::exsig::ex_double, "expected 2200 (aka fermion double) exsig");
    if (!conn.kramers_conserve()) return 0;
    // fermi phase is always negative
    return -HeisenbergFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[1], conn.m_ann[0]);
}