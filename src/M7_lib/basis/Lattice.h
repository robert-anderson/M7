//
// Created by Robert J. Anderson on 2/3/22.
//

#ifndef M7_LATTICE_H
#define M7_LATTICE_H

#include <utility>

#include <M7_lib/linalg/Dense.h>
#include <M7_lib/nd/NdFormatD.h>
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
        uint_t m_isite;
        int m_phase;
        bool operator==(const AdjElement& other) const;
    };
    typedef v_t<AdjElement> adj_row_t;

    struct Base {
        const uint_t m_nsite;
        /**
         * number of adjacent sites for each site
         */
        const uintv_t m_nadjs;
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
        const uint_t m_unique_nadj_product;
        /**
         * maximum number of adjacent sites for any one site
         */
        const uint_t m_nadj_max;

        Base(const uintv_t& nadjs);
        virtual ~Base() = default;
        virtual int phase(uint_t isite, uint_t jsite) const = 0;
        virtual void get_adj_row(uint_t isite, adj_row_t &row) const = 0;
        virtual str_t info() const = 0;

        template<typename T>
        const T* as() const {
            return dynamic_cast<const T*>(this);
        }

        operator bool() const {
            return m_nsite;
        }

    private:
        uint_t make_unique_nadj_product();
        uint_t make_nadj_max();
    };

    struct OrthoTopology {
        const NdEnumerationD m_inds;
        const v_t<int> m_bcs;
        const str_t m_info_string;

        OrthoTopology(const uintv_t &shape, const v_t<int> &bcs);

    private:
        int one_dim_phase(uint_t iind, uint_t jind, uint_t idim) const;

    public:
        uint_t isite_adj(const uintv_t &inds, uint_t idim, uint_t value) const;

        uint_t nsite() const;

        int phase(uint_t isite, uint_t jsite) const;

        void get_adj_row(uint_t isite, lattice::adj_row_t &row) const;
    };

    /**
     * use when there is no lattice structure
     */
    struct NullTopology {
        const str_t m_info_string = "null";

        uint_t isite_adj(const uintv_t &inds, uint_t idim, uint_t value) const;

        uint_t nsite() const;

        int phase(uint_t isite, uint_t jsite) const;

        void get_adj_row(uint_t isite, lattice::adj_row_t &row) const;
    };

    template<typename topo_t>
    struct Lattice : Base {
        /**
         * instance of the topology implementation
         */
        const topo_t m_topo;
    protected:
        static uintv_t make_nadjs(const topo_t& topo) {
            const auto nsite = topo.nsite();
            uintv_t out;
            out.reserve(nsite);
            adj_row_t row;
            for (uint_t isite=0; isite < nsite; ++isite) {
                topo.get_adj_row(isite, row);
                out.push_back(row.size());
            }
            return out;
        }
    public:

        Lattice(const topo_t& topo): Base(make_nadjs(topo)), m_topo(topo){}

        int phase(uint_t isite, uint_t jsite) const override {
            return m_topo.phase(isite, jsite);
        }

        void get_adj_row(uint_t isite, adj_row_t &row) const override {
            return m_topo.get_adj_row(isite, row);
        }

        str_t info() const override {
            return m_topo.m_info_string;
        }
    };

    typedef Lattice<OrthoTopology> Ortho;
    typedef Lattice<NullTopology> Null;

    /**
     * @return
     *  shared pointer to a null lattice
     */
    std::shared_ptr<Base> make();
    std::shared_ptr<Base> make(str_t topo, uintv_t site_shape, v_t<int> bcs);
    std::shared_ptr<Base> make(const conf::LatticeModel& opts);
}

#endif //M7_LATTICE_H
