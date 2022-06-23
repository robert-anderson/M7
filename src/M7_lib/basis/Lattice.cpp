//
// Created by Robert J. Anderson on 2/3/22.
//

#include "Lattice.h"
#include "M7_lib/util/SmartPtr.h"

bool lattice::AdjElement::operator==(const lattice::AdjElement &other) const {
    return m_isite==other.m_isite && m_phase==other.m_phase;
}

lattice::Base::Base(const defs::inds_t &nadjs) :
        m_nsite(nadjs.size()), m_nadjs(nadjs), m_unique_nadj_product(make_unique_nadj_product()),
        m_nadj_max(make_nadj_max()){}

size_t lattice::Base::make_unique_nadj_product() {
    size_t out = 1ul;
    for (size_t nadj: m_nadjs) {
        if (!nadj) continue;
        if ((out/nadj)*nadj != out) out*=nadj;
    }
    return out;
}

size_t lattice::Base::make_nadj_max() {
    if (!*this) return 0ul;
    REQUIRE_FALSE(m_nadjs.empty(), "number of adjacent sites vector is not set");
    return *std::max_element(m_nadjs.cbegin(), m_nadjs.cend());
}

lattice::OrthoTopology::OrthoTopology(const defs::inds_t &shape, const std::vector<int> &bcs) :
    m_inds(shape), m_bcs(bcs),
    m_info_string(log::format("orthogonal lattice with shape {} and boundary conds {}",
                              convert::to_string(m_inds.m_shape), convert::to_string(m_bcs))) {}

int lattice::OrthoTopology::one_dim_phase(size_t iind, size_t jind, size_t idim) const {
    if ((iind + 1 == jind) || (iind - 1 == jind)) return 1;
    auto bc = m_bcs[idim];
    if (!bc) return 0;
    auto max_ind = m_inds.m_shape[idim] - 1;
    if ((iind == 0 && jind == max_ind) || (jind == 0 && iind == max_ind)) return bc;
    return 0;
}

size_t lattice::OrthoTopology::isite_adj(const defs::inds_t &inds, size_t idim, size_t value) const {
    auto orig_value = inds[idim];
    auto &mutable_inds = const_cast<defs::inds_t &>(inds);
    mutable_inds[idim] = value;
    auto i = m_inds.flatten(mutable_inds);
    // leave inds_t unchanged
    mutable_inds[idim] = orig_value;
    return i;
}

size_t lattice::OrthoTopology::nsite() const {
    return m_inds.m_nelement;
}

int lattice::OrthoTopology::phase(size_t isite, size_t jsite) const {
    const auto &iinds = m_inds[isite];
    const auto &jinds = m_inds[jsite];
    int res = 0;
    const auto ndim = m_inds.m_nind;
    for (size_t idim = 0ul; idim < ndim; ++idim) {
        auto adj = one_dim_phase(iinds[idim], jinds[idim], idim);
        // only return non-zero if adjacent in a single dimension
        if (adj) {
            if (res) return 0;
            else res = adj;
        }
    }
    return res;
}

void lattice::OrthoTopology::get_adj_row(size_t isite, lattice::adj_row_t &row) const {
    row.clear();
    const auto &iinds = m_inds[isite];
    const auto ndim = m_inds.m_nind;
    for (size_t idim = 0ul; idim < ndim; ++idim) {
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

size_t lattice::NullTopology::isite_adj(const defs::inds_t &inds, size_t idim, size_t value) const {
    return ~0ul;
}

size_t lattice::NullTopology::nsite() const {
    return 0ul;
}

int lattice::NullTopology::phase(size_t isite, size_t jsite) const {
    return 0;
}

void lattice::NullTopology::get_adj_row(size_t isite, lattice::adj_row_t &row) const {
    row.clear();
}

std::shared_ptr<lattice::Base> lattice::make() {
    return smart_ptr::make_poly_shared<Base, Null>(NullTopology());
}

std::shared_ptr<lattice::Base> lattice::make(std::string topo, defs::inds_t site_shape, std::vector<int> bcs) {
    if (topo == "ortho" || topo == "orthogonal")
        return smart_ptr::make_poly_shared<Base, Ortho>(OrthoTopology(site_shape, bcs));
    return make();
}

std::shared_ptr<lattice::Base> lattice::make(const conf::LatticeModel &opts) {
    return make(opts.m_topology, opts.m_site_shape, opts.m_boundary_conds);
}
