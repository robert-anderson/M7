//
// Created by Robert J. Anderson on 2/3/22.
//

#include "Lattice.h"
#include "M7_lib/util/SmartPtr.h"

lattice::Lattice::Lattice(const sparse::dynamic::Matrix<int>& adj, uint_t nsite, str_t info_str) :
        m_nsite(nsite), m_info(std::move(info_str)), m_sparse_adj(adj), m_nadj_max(m_sparse_adj.m_max_nentry),
        m_sparse_inv(m_sparse_adj), m_lcm_le_nadj_max(integer::lcm_le(m_nadj_max)){
    REQUIRE_LE(adj.nrow(), m_nsite, "highest site index in adjacency map is OOB");
    for (uint_t isite=0ul; isite<m_nsite; ++isite) {
        if (isite >= adj.nrow() || !adj.nentry(isite))
            log::warn("lattice site {} has no adjacent sites!", isite);
    }
}

lattice::Lattice::Lattice(const lattice::Topology &topo) : Lattice(topo.make_adj(), topo.m_nsite, topo.m_info){}

lattice::OrthoTopology::OrthoTopology(const uintv_t &shape, const v_t<int> &bcs) :
        Topology(NdFormatD(shape).m_nelement,
                 log::format("orthogonal lattice with shape {} and boundary conds {}",
                             convert::to_string(m_inds.m_shape), convert::to_string(m_bcs))),
    m_inds(shape), m_bcs(bcs) {
    for (uint_t idim=0ul; idim<m_inds.m_nind; ++idim){
        REQUIRE_TRUE(m_inds.m_shape[idim], "every extent in the site shape must be non-zero");
        if (m_inds.m_shape[idim] == 1ul)
            log::warn("redundant dimension in orthogonal lattice definition");
        if (m_inds.m_shape[idim] <= 2ul)
            REQUIRE_FALSE(m_bcs[idim], "(anti-)periodic boundary conditions are not compatible with extents <= 2");
    }
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

lattice::Topology::adj_t lattice::OrthoTopology::make_adj() const {
    adj_t adj;
    adj.resize(m_nsite);
    for (uint_t isite=0ul; isite<m_inds.m_nelement; ++isite) {
        const auto &iinds = m_inds[isite];
        const auto ndim = m_inds.m_nind;
        for (uint_t idim = 0ul; idim < ndim; ++idim) {
            const auto ind = iinds[idim];
            const auto max_ind = m_inds.m_shape[idim] - 1;
            const auto bc = m_bcs[idim];
            if (!max_ind) continue;
            if (max_ind == 1ul) {
                adj.insert(isite, {isite_adj(iinds, idim, !ind), 1});
                continue;
            }
            if (ind == 0) {
                if (bc) adj.insert(isite, {isite_adj(iinds, idim, max_ind), bc});
                adj.insert(isite, {isite_adj(iinds, idim, 1ul), 1});
            } else if (ind == max_ind) {
                adj.insert(isite, {isite_adj(iinds, idim, max_ind - 1), 1});
                if (bc) adj.insert(isite, {isite_adj(iinds, idim, 0ul), bc});
            } else {
                adj.insert(isite, {isite_adj(iinds, idim, ind - 1), 1});
                adj.insert(isite, {isite_adj(iinds, idim, ind + 1), 1});
            }
        }
    }
    return adj;
}


std::shared_ptr<lattice::Lattice> lattice::make() {
    return smart_ptr::make_shared<Lattice>(NullTopology());
}

std::shared_ptr<lattice::Lattice> lattice::make(str_t topo, uintv_t site_shape, v_t<int> bcs) {
    if (topo == "ortho" || topo == "orthogonal")
        return smart_ptr::make_shared<Lattice>(OrthoTopology(site_shape, bcs));
    return make();
}

std::shared_ptr<lattice::Lattice> lattice::make(const conf::LatticeModel &opts) {
    return make(opts.m_topology, opts.m_site_shape, opts.m_boundary_conds);
}
