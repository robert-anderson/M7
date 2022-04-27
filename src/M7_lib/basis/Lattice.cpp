//
// Created by anderson on 2/3/22.
//

#include "Lattice.h"

Lattice::Topology Lattice::topology(const std::string &str) {
    if (str=="ortho") return Ortho;
    if (str=="orthogonal") return Ortho;
    ABORT("invalid name given for lattice topology");
    return NullTopology;
}

std::string Lattice::topo_string(Lattice::Topology topo) {
    switch (topo) {
        case Ortho: return "orthogonal";
        case NullTopology: return "null";
    }
    return {};
}

const size_t &Lattice::nsite() const {
    return m_dense.nrow();
}

std::string Lattice::info() const {
    return log::format("{} lattice with site shape: {}, boundary conds {}",
                       topo_string(m_spec.m_topo), utils::to_string(m_spec.m_format.m_shape),
                       utils::to_string(m_spec.m_bcs));
}

void Lattice::add(size_t irow, const defs::inds &icols, const std::vector<int> &coeffs) {
    REQUIRE_TRUE(m_sparse.empty(irow), "row already added");
    REQUIRE_EQ(icols.size(), coeffs.size(), "unequal lengths of column indices and values");
    auto ncol = icols.size();
    if ((m_unique_nconn_product/ncol)*ncol != m_unique_nconn_product) m_unique_nconn_product*=ncol;
    for (size_t iicol=0ul; iicol<ncol; ++iicol) {
        auto icol = icols[iicol];
        auto coeff = coeffs[iicol];
        m_sparse.add(irow, {icol, coeff});
        m_dense(irow, icol) = coeff;
    }
}

Lattice::Lattice(Lattice::Spec spec) :
        m_dense(spec.m_format.m_nelement), m_spec(std::move(spec)) {
    m_sparse.resize(nsite());
}

Lattice::Spec::Spec(Lattice::Topology topo, const defs::inds &shape, std::vector<int> bcs) :
        m_topo(topo), m_format(shape), m_bcs(std::move(bcs)){
    REQUIRE_EQ(m_format.m_nind, m_bcs.size(),
               "site shape and boundary conds should be the same length");
}

Lattice::Spec::Spec(const std::string &topo, const defs::inds &shape, std::vector<int> bcs) :
        Spec(topology(topo), shape, std::move(bcs)){}

size_t OrthoLattice::get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const {
    DEBUG_ASSERT_LT(value, m_spec.m_format.m_shape[idim], "site inds value OOB");
    auto orig_value = site_inds[idim];
    auto &inds = const_cast<defs::inds &>(site_inds);
    inds[idim] = value;
    auto i = m_spec.m_format.flatten(inds);
    // leave inds unchanged
    inds[idim] = orig_value;
    return i;
}

std::pair<size_t, int> OrthoLattice::get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const {
    auto &format = m_spec.m_format;
    auto &bcs = m_spec.m_bcs;
    size_t dim_ind = ~0ul;
    int sign = 0;
    if (!inc && site_inds[idim] == 0) {
        // lower boundary
        if (bcs[idim]) {
            dim_ind = format.m_shape[idim] - 1;
            sign = bcs[idim];
        }
    } else if (inc && (site_inds[idim] + 1 == format.m_shape[idim])) {
        // upper boundary
        if (bcs[idim]) {
            dim_ind = 0ul;
            sign = bcs[idim];
        }
    } else {
        // not at a boundary
        dim_ind = site_inds[idim] + (inc ? 1 : -1);
        sign = 1;
    }
    if (dim_ind == ~0ul) return {~0ul, 0};
    return {get_coord_index(site_inds, idim, dim_ind), sign};
}

OrthoLattice::OrthoLattice(Lattice::Spec spec) : Lattice(spec){
    /*
     * set up loop for orthogonally-coordinated lattice
     */
    foreach::rtnd::Unrestricted loop(m_spec.m_format.m_shape);
    defs::inds cols;
    std::vector<int> signs;
    auto fn = [&]() {
        auto &inds = loop.m_inds;
        cols.clear();
        signs.clear();
        for (size_t idim = 0ul; idim < inds.size(); ++idim) {
            for (bool inc : {false, true}){
                auto pair = get_coordination(inds, idim, inc);
                if (pair.first != ~0ul) {
                    cols.push_back(pair.first);
                    signs.push_back(pair.second);
                }
            }
        }
        auto irow = m_spec.m_format.flatten(inds);
        add(irow, cols, signs);
    };
    loop(fn);
}

Lattice lattice::make(const Lattice::Spec &spec) {
    switch (spec.m_topo) {
        case Lattice::Ortho:
            return OrthoLattice(spec);
        case Lattice::NullTopology:
            ABORT("invalid lattice topology given");
    }
    return {Lattice::Spec(Lattice::NullTopology, {}, {})};
}

Lattice lattice::make(const conf::LatticeModel &opts) {
    return make({opts.m_topology, opts.m_site_shape, opts.m_boundary_conds});
}
