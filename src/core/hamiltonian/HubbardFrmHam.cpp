//
// Created by anderson on 12/8/21.
//

#include "HubbardFrmHam.h"


bool HubbardFrmHam::sign_problem() const {
    // can only be SPF in 1D
    if (m_format.m_nind!=1) return false;
    // open boundary conditions is always SPF
    auto bc = m_bcs[0];
    if (!bc) return true;
    auto odd_nalpha = ci_utils::nalpha(m_nelec, m_ms2_restrict) & 1ul;
    auto odd_nbeta = ci_utils::nbeta(m_nelec, m_ms2_restrict) & 1ul;
    // we have (anti-)periodic BCs: if alpha and beta oddness is different, there is a sign problem
    if (odd_nalpha != odd_nbeta) return false;
    // if nalpha is odd, then boundary excitation does not pick up a factor of -1 from the fermi phase => PBC
    // if nalpha is even, then boundary excitation does pick up a factor of -1 from the fermi phase => APBC
    return odd_nalpha == (bc > 0);
}

HubbardFrmHam::HubbardFrmHam(defs::ham_t u, Lattice lattice, int ms2_restrict, int charge) :
        FrmHam(lattice.nsite()-charge, lattice.nsite(), ms2_restrict),
        m_u(u), m_lattice(std::move(lattice)), m_format(m_lattice.m_spec.m_format),
        m_bcs(m_lattice.m_spec.m_bcs), m_spf(sign_problem()){
    m_contribs_1100.set_nonzero(exsig_utils::ex_single);
    m_contribs_2200.set_nonzero(0);
    log::info("Hubbard Hamiltonian initialized with U={}; {}", m_u, m_lattice.info());
    log::info("This model {} sign problem-free", m_spf ? "is" : "is not");
}

HubbardFrmHam::HubbardFrmHam(const fciqmc_config::FermionHamiltonian &opts) :
        HubbardFrmHam(opts.m_hubbard.m_repulsion,
                      lattice::make(Lattice::Spec(
                              opts.m_hubbard.m_topology, opts.m_hubbard.m_site_shape,
                              opts.m_hubbard.m_boundary_conds)), opts.m_ms2_restrict, opts.m_charge){}

defs::ham_t HubbardFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t h = 0.0;
    for (size_t isite = 0ul; isite < m_nsite; ++isite)
        if (onv.get({0, isite}) && onv.get({1, isite})) h += m_u;
    return h;
}

defs::ham_t HubbardFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.size(), 2ul, "incorrect connection exsig");
    auto isite = onv.isite(conn.m_ann[0]);
    auto jsite = onv.isite(conn.m_cre[0]);
    int t_mat_element = m_lattice.m_dense(isite, jsite);
    if (!t_mat_element) return 0.0;
    // don't need to compute fermi phase if the model meets the SPF conditions
    if (m_spf) return -1;
    return conn.phase(onv) ? -t_mat_element : t_mat_element;
}

defs::ham_t HubbardFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    return 0;
}

void HubbardFrmHam::log_data() const {
    FrmHam::log_data();
}

defs::ham_t HubbardFrmHam::get_coeff_1100(const size_t &i, const size_t &j) const {
    return m_lattice.m_dense(i, j);
}

defs::ham_t HubbardFrmHam::get_coeff_2200(const size_t &i, const size_t &j,
                                               const size_t &k, const size_t &l) const {
    return 0.0;
}
