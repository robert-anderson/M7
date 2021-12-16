//
// Created by anderson on 12/8/21.
//

#include "HubbardFrmHam.h"

size_t HubbardFrmHam::get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const {
    DEBUG_ASSERT_LT(value, m_format.m_shape[idim], "site inds value OOB");
    auto orig_value = site_inds[idim];
    auto &inds = const_cast<defs::inds &>(site_inds);
    inds[idim] = value;
    auto i = m_format.flatten(inds);
    // leave inds unchanged
    inds[idim] = orig_value;
    return i;
}

std::pair<size_t, int> HubbardFrmHam::get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const {
    size_t dim_ind = ~0ul;
    int t_element = 0;
    if (!inc && site_inds[idim] == 0) {
        // lower boundary
        if (m_bcs[idim]) {
            dim_ind = m_format.m_shape[idim] - 1;
            t_element = -m_bcs[idim];
        }
    } else if (inc && (site_inds[idim] + 1 == m_format.m_shape[idim])) {
        // upper boundary
        if (m_bcs[idim]) {
            dim_ind = 0ul;
            t_element = -m_bcs[idim];
        }
    } else {
        // not at a boundary
        dim_ind = site_inds[idim] + (inc ? 1 : -1);
        t_element = -1;
    }
    if (dim_ind == ~0ul) return {~0ul, 0};
    return {get_coord_index(site_inds, idim, dim_ind), t_element};
}

size_t HubbardFrmHam::nsite(const defs::inds &site_shape) {
    return NdFormatD(site_shape).m_nelement;
}

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



HubbardFrmHam::HubbardFrmHam(const defs::inds& site_shape, const std::vector<int>& bcs, defs::ham_t u,
        int ms2_restrict, int charge) :
        FrmHam(nsite(site_shape)-charge, nsite(site_shape), ms2_restrict),
        m_format(site_shape), m_bcs(bcs), m_u(u), m_t_mat_dense(m_nsite), m_spf(sign_problem()){

    m_contribs_1100.set_nonzero(exsig_utils::ex_single);
    m_contribs_2200.set_nonzero(0);
    REQUIRE_EQ(site_shape.size(), bcs.size(), "site shape and boundary conds should be the same length");

    m_t_mat_sparse.resize(m_nsite);
    m_t_mat_dense.zero();

    foreach::rtnd::Unrestricted loop(m_format.m_shape);
    std::set<size_t> nconns;
    auto fn = [&]() {
        auto &inds = loop.m_inds;
        auto irow = m_format.flatten(inds);
        size_t nconn = 0ul;
        for (size_t idim = 0ul; idim < inds.size(); ++idim) {
            for (bool inc : {false, true}){
                auto pair = get_coordination(inds, idim, inc);
                if (pair.first != ~0ul) {
                    m_t_mat_sparse.add(irow, pair);
                    m_t_mat_dense(irow, pair.first) = pair.second;
                    ++nconn;
                }
            }
        }
        nconns.insert(nconn);
    };
    loop(fn);
    m_unique_nconn_product = 1ul;
    for (const auto& unique_nconn: nconns) m_unique_nconn_product*=unique_nconn;

    log::info("Hubbard Hamiltonian initialized with U={}, site shape={}, boundary conds={}",
              m_u, utils::to_string(m_format.m_shape), utils::to_string(m_bcs));
    log::info("This model {} sign problem-free", m_spf ? "is" : "is not");
}

HubbardFrmHam::HubbardFrmHam(const fciqmc_config::FermionHamiltonian &opts) :
        HubbardFrmHam(opts.m_hubbard.m_site_shape, opts.m_hubbard.m_boundary_conds,
                           opts.m_hubbard.m_repulsion, opts.m_ms2_restrict, opts.m_charge){}

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
    int t_mat_element = m_t_mat_dense(isite, jsite);
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
    return m_t_mat_dense(i, j);
}

defs::ham_t HubbardFrmHam::get_coeff_2200(const size_t &i, const size_t &j,
                                               const size_t &k, const size_t &l) const {
    return 0.0;
}
