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
            logging::warn("lattice site {} has no adjacent sites!", isite);
    }
}

lattice::Lattice::Lattice(const lattice::Topology &topo) : Lattice(topo.make_adj(), topo.m_nsite, topo.m_info){}

lattice::OrthoTopology::OrthoTopology(const uintv_t &shape, const v_t<int> &bcs) :
        Topology(NdFormatD(shape).m_nelement,
                 logging::format("orthogonal lattice with shape {} and boundary conds {}",
                             convert::to_string(shape), convert::to_string(bcs))),
    m_inds(shape), m_bcs(bcs) {
    for (uint_t idim=0ul; idim<m_inds.m_nind; ++idim){
        REQUIRE_TRUE(m_inds.m_shape[idim], "every extent in the site shape must be non-zero");
        if (m_inds.m_shape[idim] == 1ul)
            logging::warn("redundant dimension in orthogonal lattice definition");
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
    /*
     * number of dimensions in the lattice
     */
    const auto ndim = m_inds.m_nind;
    for (uint_t isite=0ul; isite<m_inds.m_nelement; ++isite) {
        /*
         * ndim-dimensional coordinates of isite
         */
        const auto &iinds = m_inds[isite];
        for (uint_t idim = 0ul; idim < ndim; ++idim) {
            /*
             * inserts an (index, phase) pair into the current isite's adjacency matrix row
             */
            auto insert = [&](uint_t jsite, int phase) {
                DEBUG_ASSERT_TRUE(phase, "shouldn't add lattice unconnected pairs");
                sparse::MatrixElement<int> adj_elem(isite_adj(iinds, idim, jsite), phase);
                adj.insert(isite, adj_elem);
            };
            /*
             * number of sites along this dimension of the lattice
             */
            const auto ind = iinds[idim];
            /*
             * maximum value taken by the indices along this dimension
             */
            const auto max_ind = m_inds.m_shape[idim] - 1;
            // go to next dimension if this one only contains a single site
            if (!max_ind) continue;
            if (max_ind == 1ul) {
                // if this dimension only has two sites, then they are set to adjacent and BCs are not considered
                insert(!ind, 1);
            }
            else {
                const auto bc = m_bcs[idim];
                if (ind == 0) {
                    // lower boundary in dimension idim
                    if (bc) insert(max_ind, bc);
                    insert(1ul, 1);
                } else if (ind == max_ind) {
                    // upper boundary in dimension idim
                    insert(max_ind - 1, 1);
                    if (bc) insert(0ul, bc);
                } else {
                    // not a boundary site
                    insert(ind - 1, 1);
                    insert(ind + 1, 1);
                }
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
