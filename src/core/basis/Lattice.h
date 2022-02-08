//
// Created by anderson on 2/3/22.
//

#ifndef M7_LATTICE_H
#define M7_LATTICE_H

#include <src/core/linalg/Dense.h>
#include <src/core/nd/NdFormatD.h>

#include <utility>
#include <src/core/util/Foreach.h>
#include <src/core/config/Hamiltonian.h>
#include "src/core/linalg/Sparse.h"

/**
 * In lattice models and potentially other systems of interest, the sums over sites, modes, or spin orbitals in the
 * Hamiltonian definition are often strictly constrained, e.g. to nearest neighbors only. Such constraints are dealt
 * with in this class in both sparse and dense representations.
 *
 * This class is not limited to dealing with nearest neighbors - any constrained sum over a pair of degrees of freedom
 * can be efficiently represented here
 */
struct Lattice {
    /**
     * sparse map enabling lookup of all coordinated sites given a row site
     */
    sparse::Matrix<int> m_sparse;
    /**
     * dense map of coordinated sites allowing lookup of the 1-body H matrix element given the flat indices of two site
     * index vectors
     */
    dense::SquareMatrix<int> m_dense;

    enum Topology {Ortho, NullTopology};

private:
    static Topology topology(const std::string& str);

    static std::string topo_string(Topology topo);

public:
    struct Spec {
        Topology m_topo;
        NdFormatD m_format;
        std::vector<int> m_bcs;
        Spec(Topology topo, const defs::inds& shape, std::vector<int> bcs);
        Spec(const std::string& topo, const defs::inds& shape, std::vector<int> bcs);
    };
    /**
     * complete specification for all lattice subtypes
     */
    const Spec m_spec;
    /**
     * to help with excitation generation. The number of connections possible from a site i cannot be known before i is
     * selected, and it is efficient to select i and a connected j in a single random number draw. If the random number
     * R is chosen is in the range [0, nsite*m_unique_nconn_product) then R/m_unique_nconn_product gives i and also
     * therefore the number of connections site i actually has (nconn_i) can be looked up via this class. But since the
     * range is an integral multiple of nconn_i, the modular remainder R%nconn_i is a suitable selecting index for the
     * connected site.
     *
     * when a set of connections is reached for which is not a divisor of this member, then it is multiplied by the new
     * number of connections. Thus the ultimate value of this integer is divisible by all nconns
     */
    size_t m_unique_nconn_product = 1ul;

    const size_t& nsite() const;

    std::string info() const;

    void add(size_t irow, const defs::inds& icols, const std::vector<int>& coeffs);

    Lattice(size_t nsite, Spec spec);

};

/**
 * 1D chains, 2D square lattices and higher-dimensional analogs thereof
 */
class OrthoLattice : public Lattice {

    size_t get_coord_index(const defs::inds &site_inds, size_t idim, size_t value) const;

    std::pair<size_t, int> get_coordination(const defs::inds &site_inds, size_t idim, bool inc) const;

public:
    OrthoLattice(Spec spec);
};

namespace lattice {

    Lattice make(const Lattice::Spec &spec);

    Lattice make(const fciqmc_config::LatticeModel& opts);
}
#endif //M7_LATTICE_H
