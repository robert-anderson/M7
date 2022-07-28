//
// Created by Robert J. Anderson on 2/3/22.
//

#ifndef M7_LATTICE_H
#define M7_LATTICE_H

#include <utility>

#include <M7_lib/linalg/Dense.h>
#include <M7_lib/nd/NdFormatD.h>
#include <M7_lib/conf/HamiltonianConf.h>
#include "M7_lib/linalg/sparse/Inverse.h"

/**
 * In lattice models and potentially other systems of interest, the sums over sites, modes, or spin orbitals in the
 * Hamiltonian definition are often strictly constrained, e.g. to nearest neighbors only. Such constraints are dealt
 * with in this class in both sparse and dense representations.
 *
 * This class is not limited to dealing with nearest neighbors - any constrained sum over a pair of degrees of freedom
 * can be efficiently represented here
 */
namespace lattice {

    struct Topology {
        const uint_t m_nsite;
        const str_t m_info;
        Topology(uint_t nsite, str_t info): m_nsite(nsite), m_info(std::move(info)){}
        typedef sparse::dynamic::Matrix<int> adj_t;
        virtual adj_t make_adj() const = 0;
        virtual ~Topology() = default;
    };

    struct SubLattice {
        /**
         * number of sites in the basis
         */
        const uint_t m_nsite;
        /**
         * information string
         */
         const str_t m_info;
        /**
         * sparse map of adjacent sites (access to connected sites and phases given a site index)
         */
        typedef sparse::fixed::Matrix<int> sparse_adj_t;
        const sparse_adj_t m_sparse_adj;
        typedef sparse_adj_t::elem_t adj_t;

        /**
         * maximum number of adjacent sites for any one site
         */
        const uint_t m_nadj_max;
        /**
         * inverse sparse map (access to adjacency + phase given two site indices)
         */
        typedef sparse::inverse::Matrix<int> sparse_inv_t;
        const sparse_inv_t m_sparse_inv;
        /**
         * to help with excitation generation. The number of connections possible from a site i cannot be known before i is
         * selected, and it is efficient to select i and a connected j in a single random number draw. If the random number
         * R is chosen is in the range [0, nsite*m_unique_nconn_product) then R/m_unique_nconn_product gives i and also
         * therefore the number of connections site i actually has (nconn_i) can be looked up via this class. But since the
         * range is an integral multiple of nconn_i, the modular remainder R%nconn_i is a suitable selecting index for the
         * connected site.
         *
         * therefore, the least common multiple of all whole numbers less than or equal to the maximum number of
         * adjacent sites is stored to facilitate this efficient drawing
         */
        const uint_t m_lcm_le_nadj_max;

        operator bool() const {
            return m_nsite;
        }

    protected:
        SubLattice(const sparse::dynamic::Matrix<int>& adj, uint_t nsite, str_t info_str);

        /**
         * a next-nearest neighbor of a lattice site i is a site j which is unconnected to i but has the maximum number
         * of mutually adjacent sites with i
         * @return
         */
        SubLattice make_next_nearest() const;

        explicit SubLattice(const Topology& topo);
    };

    struct Lattice : SubLattice {
        /**
         * lattice of next-nearest neighbors
         */
        const SubLattice m_next_nearest;

        explicit Lattice(const Topology& topo): SubLattice(topo), m_next_nearest(make_next_nearest()){}
    };

    struct OrthoTopology : Topology {
        const NdEnumerationD m_inds;
        const v_t<int> m_bcs;

        OrthoTopology(const uintv_t &shape, const v_t<int> &bcs);

    private:
        uint_t isite_adj(const uintv_t &inds, uint_t idim, uint_t value) const;
    public:

        adj_t make_adj() const override;
    };

    /**
     * use when there is no lattice structure
     */
    struct NullTopology : Topology {

        NullTopology() : Topology(0ul, "null") {}

        adj_t make_adj() const override {
            return {};
        }
    };


    /**
     * @return
     *  shared pointer to a null lattice
     */
    std::shared_ptr<Lattice> make();
    std::shared_ptr<Lattice> make(str_t topo, uintv_t site_shape, v_t<int> bcs);
    std::shared_ptr<Lattice> make(const conf::LatticeModel& opts);
}

#endif //M7_LATTICE_H
