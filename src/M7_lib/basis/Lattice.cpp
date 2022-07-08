//
// Created by Robert J. Anderson on 2/3/22.
//

#include "Lattice.h"
#include "M7_lib/util/SmartPtr.h"

bool lattice::AdjElement::operator==(const lattice::AdjElement &other) const {
    return m_isite==other.m_isite && m_phase==other.m_phase;
}

lattice::Base::Base(const uintv_t &nadjs) :
        m_nsite(nadjs.size()), m_nadjs(nadjs), m_nadj_max(make_nadj_max()),
        m_lcm_le_nadj_max(integer::lcm_le(m_nadj_max)){}
        
uint_t lattice::Base::make_nadj_max() {
    if (!*this) return 0ul;
    REQUIRE_FALSE(m_nadjs.empty(), "number of adjacent sites vector is not set");
    return *std::max_element(m_nadjs.cbegin(), m_nadjs.cend());
}

lattice::OrthoTopology::OrthoTopology(const uintv_t &shape, const v_t<int> &bcs) :
    m_inds(shape), m_bcs(bcs),
    m_info_string(log::format("orthogonal lattice with shape {} and boundary conds {}",
                              convert::to_string(m_inds.m_shape), convert::to_string(m_bcs))) {}

int lattice::OrthoTopology::one_dim_phase(uint_t iind, uint_t jind, uint_t idim) const {
    if ((iind + 1 == jind) || (iind - 1 == jind)) return 1;
    auto bc = m_bcs[idim];
    if (!bc) return 0;
    auto max_ind = m_inds.m_shape[idim] - 1;
    if ((iind == 0 && jind == max_ind) || (jind == 0 && iind == max_ind)) return bc;
    return 0;
}

uint_t lattice::OrthoTopology::isite_adj(const uintv_t &inds, uint_t idim, uint_t value) const {
    auto orig_value = inds[idim];
    auto &mutable_inds = const_cast<uintv_t &>(inds);
    mutable_inds[idim] = value;
    auto i = m_inds.flatten(mutable_inds);
    // leave uintv_t unchanged
    mutable_inds[idim] = orig_value;
    return i;
}

uint_t lattice::OrthoTopology::nsite() const {
    return m_inds.m_nelement;
}

int lattice::OrthoTopology::phase(uint_t isite, uint_t jsite) const {
    const auto &iinds = m_inds[isite];
    const auto &jinds = m_inds[jsite];
    int res = 0;
    const auto ndim = m_inds.m_nind;
    for (uint_t idim = 0ul; idim < ndim; ++idim) {
        auto adj = one_dim_phase(iinds[idim], jinds[idim], idim);
        // only return non-zero if adjacent in a single dimension
        if (adj) {
            if (res) return 0;
            else res = adj;
        }
    }
    return res;
}

void lattice::OrthoTopology::get_adj_row(uint_t isite, lattice::adj_row_t &row) const {
    row.clear();
    const auto &iinds = m_inds[isite];
    const auto ndim = m_inds.m_nind;
    for (uint_t idim = 0ul; idim < ndim; ++idim) {
        const auto ind = iinds[idim];
        const auto max_ind = m_inds.m_shape[idim] - 1;
        const auto bc = m_bcs[idim];
        if (!max_ind) continue;
        if (max_ind == 1ul) {
            row.push_back({isite_adj(iinds, idim, !ind), 1});
            continue;
        }
        if (ind == 0) {
            if (bc) row.push_back({isite_adj(iinds, idim, max_ind), bc});
            row.push_back({isite_adj(iinds, idim, 1ul), 1});
        }
        else if (ind == max_ind) {
            row.push_back({isite_adj(iinds, idim, max_ind - 1), 1});
            if (bc) row.push_back({isite_adj(iinds, idim, 0ul), bc});
        }
        else {
            row.push_back({isite_adj(iinds, idim, ind - 1), 1});
            row.push_back({isite_adj(iinds, idim, ind + 1), 1});
        }
    }
}

uint_t lattice::NullTopology::isite_adj(const uintv_t &/*inds*/, uint_t /*idim*/, uint_t /*value*/) const {
    return ~0ul;
}

uint_t lattice::NullTopology::nsite() const {
    return 0ul;
}

int lattice::NullTopology::phase(uint_t /*isite*/, uint_t /*jsite*/) const {
    return 0;
}

void lattice::NullTopology::get_adj_row(uint_t /*isite*/, lattice::adj_row_t &row) const {
    row.clear();
}

std::shared_ptr<lattice::Base> lattice::make() {
    return smart_ptr::make_poly_shared<Base, Null>(NullTopology());
}

std::shared_ptr<lattice::Base> lattice::make(str_t topo, uintv_t site_shape, v_t<int> bcs) {
    if (topo == "ortho" || topo == "orthogonal")
        return smart_ptr::make_poly_shared<Base, Ortho>(OrthoTopology(site_shape, bcs));
    return make();
}

std::shared_ptr<lattice::Base> lattice::make(const conf::LatticeModel &opts) {
    return make(opts.m_topology, opts.m_site_shape, opts.m_boundary_conds);
}
