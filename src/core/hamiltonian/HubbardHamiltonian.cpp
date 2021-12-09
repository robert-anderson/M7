//
// Created by anderson on 12/8/21.
//

#include "HubbardHamiltonian.h"

size_t HubbardHamiltonian::get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const {
    auto orig_value = site_inds[idim];
    auto &inds = const_cast<defs::inds &>(site_inds);
    inds[idim] = value;
    auto i = m_format.flatten(inds);
    // leave inds unchanged
    inds[idim] = orig_value;
    return i;
}

std::pair<size_t, int> HubbardHamiltonian::get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const {
    size_t dim_ind = ~0ul;
    int t_element = 0;
    if (site_inds[idim] == 0 && !inc) {
        // lower boundary
        if (m_bcs[idim]) {
            dim_ind = m_format.m_shape[idim] - 1;
            t_element = -m_bcs[idim];
        }
    } else if (site_inds[idim] + 1 == m_format.m_shape[idim] && inc) {
        // upper boundary
        if (!m_bcs[idim]) {
            dim_ind = 0ul;
            t_element = -m_bcs[idim];
        }
    } else {
        // not at a boundary
        dim_ind = inc ? dim_ind + 1 : dim_ind - 1;
        t_element = -1;
    }
    if (dim_ind == ~0ul) return {~0ul, 0};
    return {get_coord_index(site_inds, idim, dim_ind), t_element};
}

size_t HubbardHamiltonian::nsite(const defs::inds &site_shape) {
    return NdFormatD(site_shape).m_nelement;
}

HubbardHamiltonian::HubbardHamiltonian(const defs::inds& site_shape, std::vector<int> bcs, defs::ham_t u,
        int ms2_restrict, int charge) :
        FermionHamiltonian(nsite(site_shape)-charge, nsite(site_shape), ms2_restrict),
        m_format(site_shape), m_u(u), m_t_mat_dense(m_nsite) {
    m_t_mat_sparse.resize(m_nsite);
    m_t_mat_dense.zero();
    foreach::rtnd::Unrestricted loop(m_format.m_shape);
    auto fn = [&]() {
        auto &inds = loop.m_inds;
        auto irow = m_format.flatten(inds);
        for (size_t idim = 0ul; idim < inds.size(); ++idim) {
            auto pair = get_coordination(inds, idim, false);
            if (pair.first != ~0ul) {
                m_t_mat_sparse.add(irow, pair);
                m_t_mat_dense(irow, pair.first) = pair.second;
            }
            pair = get_coordination(inds, idim, true);
            if (pair.first != ~0ul) {
                m_t_mat_sparse.add(irow, pair);
                m_t_mat_dense(irow, pair.first) = pair.second;
            }
        }
    };
    loop(fn);
}

defs::ham_t HubbardHamiltonian::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t h = 0.0;
    for (size_t isite = 0ul; isite < m_nsite; ++isite)
        if (onv.get({0, isite}) && onv.get({1, isite})) h += m_u;
    return h;
}

defs::ham_t HubbardHamiltonian::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.size(), 2ul, "incorrect connection exsig");
    int t_mat_element = m_t_mat_dense(conn.m_ann[0], conn.m_cre[0]);
    if (!t_mat_element) return 0.0;
    return conn.phase(onv) ? -t_mat_element : t_mat_element;
}

defs::ham_t HubbardHamiltonian::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    return 0;
}

void HubbardHamiltonian::log_data() const {
    FermionHamiltonian::log_data();
}
