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

Lattice::Lattice(size_t nsite, Lattice::Spec spec) :
        m_dense(nsite), m_spec(std::move(spec)) {
    m_sparse.resize(nsite);
}

Lattice::Spec::Spec(Lattice::Topology topo, const defs::inds &shape, std::vector<int> bcs) :
        m_topo(topo), m_format(shape), m_bcs(std::move(bcs)){
    REQUIRE_EQ(m_format.m_nind, m_bcs.size(),
               "site shape and boundary conds should be the same length");
}

Lattice::Spec::Spec(const std::string &topo, const defs::inds &shape, std::vector<int> bcs) :
        Spec(topology(topo), shape, std::move(bcs)){}
