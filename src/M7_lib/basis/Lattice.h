//
// Created by Robert J. Anderson on 2/3/22.
//

#ifndef M7_LATTICE_H
#define M7_LATTICE_H

#include <utility>

#include <M7_lib/linalg/Dense.h>
#include <M7_lib/nd/NdFormatD.h>
#include <M7_lib/foreach/Foreach.h>
#include <M7_lib/conf/HamiltonianConf.h>
#include <M7_lib/linalg/Sparse.h>

/**
 * In lattice models and potentially other systems of interest, the sums over sites, modes, or spin orbitals in the
 * Hamiltonian definition are often strictly constrained, e.g. to nearest neighbors only. Such constraints are dealt
 * with in this class in both sparse and dense representations.
 *
 * This class is not limited to dealing with nearest neighbors - any constrained sum over a pair of degrees of freedom
 * can be efficiently represented here
 */
namespace lattice {

    struct AdjElement {
        const size_t m_isite;
        const int m_phase;
    };
    typedef std::vector<AdjElement> adj_row_t;

    struct Base {
        const size_t m_nsite;
        /**
         * number of adjacent sites for each site
         */
        const defs::inds m_nadjs;
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
        const size_t m_unique_nadj_product;
        /**
         * maximum number of adjacent sites for any one site
         */
        const size_t m_nadj_max;

        Base(const defs::inds& nadjs);
        virtual ~Base() = default;
        virtual int phase(size_t isite, size_t jsite) const = 0;
        virtual void get_adj_row(size_t isite, adj_row_t &row) const = 0;
        virtual std::string info() const = 0;

        template<typename T>
        const T* as() const {
            return dynamic_cast<const T*>(this);
        }

    private:
        size_t make_unique_nadj_product();
        size_t make_nadj_max();
    };


    template<typename topo_t>
    struct Lattice : Base {
        /**
         * instance of the topology implementation
         */
        const topo_t m_topo;
    protected:
        static defs::inds make_nadjs(const topo_t& topo) {
            const auto nsite = topo.nsite();
            defs::inds out;
            out.reserve(nsite);
            adj_row_t row;
            for (size_t isite=0; isite < nsite; ++isite) {
                topo.get_adj_row(isite, row);
                out.push_back(row.size());
            }
            return out;
        }
    public:

        Lattice(const topo_t& topo): Base(make_nadjs(topo)), m_topo(topo){}

        int phase(size_t isite, size_t jsite) const override {
            return m_topo.phase(isite, jsite);
        }

        void get_adj_row(size_t isite, adj_row_t &row) const override {
            return m_topo.get_adj_row(isite, row);
        }

        std::string info() const override {
            return m_topo.m_info_string;
        }
    };

    struct OrthoTopology {
        const NdEnumerationD m_inds;
        const std::vector<int> m_bcs;
        const std::string m_info_string;

        OrthoTopology(const defs::inds &shape, const std::vector<int> &bcs);

        int one_dim_phase(size_t iind, size_t jind, size_t idim) const;

        size_t isite_adj(const defs::inds &inds, size_t idim, size_t value) const;

        size_t nsite() const;

        int phase(size_t isite, size_t jsite) const;

        void get_adj_row(size_t isite, lattice::adj_row_t &row) const;
    };

    typedef Lattice<OrthoTopology> Ortho;

    std::shared_ptr<Base> make(std::string topo, defs::inds site_shape, std::vector<int> bcs);
    std::shared_ptr<Base> make(const conf::LatticeModel& opts);
}

#endif //M7_LATTICE_H
