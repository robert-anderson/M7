//
// Created by rja on 24/05/2021.
//

#ifndef M7_HAMILTONIANPARTS_H
#define M7_HAMILTONIANPARTS_H

#include <src/core/table/BufferedFields.h>
#include <src/core/io/Options.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "Hamiltonian.h"

#if 0
namespace ham_parts {

    struct Fermion {
        const size_t m_nelec;
        const size_t m_nsite;
        const bool m_spin_conserving_1e, m_spin_conserving_2e;
        const bool m_complex_valued;
        const size_t m_int_2e_rank;

    public:
        defs::ham_t m_int_0;
        typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
        typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
        ints1_t m_int_1;
        ints2_t m_int_2;

        /**
         * stores whether connected elements are nonzero based on contrib case
         */
        std::array<bool, ham_data::contrib_count> m_nonzero_contribs;
        /**
         * are (1, 1)-rank contribs due to (1, 1)-excitations nonzero only for (1D) nearest neighbors?
         * assume this is the case unless counterexample found in loop over nonzero elements
         */
        bool m_nn_only_1111 = true;
        /**
         * same as above, but with (1D) periodic boundary conditions
         */
        bool m_nnp_only_1111 = true;
        /**
         * are (2, 2)-rank contribs due to (0, 0)-excitations nonzero only when all operator indices refer to the same site?
         * assume this is the case unless counterexample found in loop over nonzero elements
         */
        bool m_on_site_only_0022 = true;

    public:

        /**
         * @param iorb
         *  orbital index (site if integrals spin resolved, else spin oribtal index)
         * @return
         *  site index
         */
        size_t iorb_to_isite(const size_t &iorb) const {
            return iorb < m_nsite ? iorb : iorb + m_nsite;
        }

        /**
         * @param iorb
         * @param jorb
         * @param periodic
         *  if true, lowest and highest site indices count as neighbors.
         * @return
         * true if orbitals are nearest neighbors
         */
        bool nearest_neighbors(const size_t &iorb, const size_t &jorb, bool periodic) const {
            auto isite = iorb_to_isite(iorb);
            auto jsite = iorb_to_isite(jorb);
            if (isite + 1 == jsite || jsite + 1 == isite) return true;
            if (periodic) {
                return (isite == 0 && jsite == m_nsite-1) || (isite == m_nsite-1 && jsite == 0);
            }
            return false;
        }

        bool on_site(const size_t &iorb, const size_t &jorb, const size_t &korb, const size_t &lorb) const {
            auto isite = iorb_to_isite(iorb);
            auto jsite = iorb_to_isite(jorb);
            auto ksite = iorb_to_isite(korb);
            auto lsite = iorb_to_isite(lorb);
            return isite == jsite && jsite == ksite && ksite == lsite;
        }

        Fermion(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                           bool complex_valued, bool spin_resolved, size_t int_2e_rank);

        Fermion(const FcidumpFileReader &file_reader);

        Fermion(std::string fname, bool spin_major);

        Fermion(const fciqmc_config::Hamiltonian &opts) : Fermion(opts.m_fcidump.m_path, opts.m_fcidump.m_spin_major) {}

        defs::ham_comp_t get_energy(const fields::Onv<0> &fonv) const;

        defs::ham_t get_element_0(const defs::inds &occs, const size_t &nocc) const;

        defs::ham_t get_element_0(const OccupiedOrbitals &occs) const;

        defs::ham_t get_element_0(const fields::Onv<0> &fonv) const;

        defs::ham_t get_element_0(const conn::Antisym<0> &connection) const;

        defs::ham_t get_element_1(const conn::Antisym<0> &connection) const;

        defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;

        defs::ham_t get_element_2(const conn::Basic<0> &connection) const;

        defs::ham_t get_element_2(const conn::Antisym<0> &connection) const;

        defs::ham_t get_element(const conn::Antisym<0> &connection) const;

        defs::ham_t get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const;

        bool is_hubbard() const {
            return m_on_site_only_0022 && m_nn_only_1111;
        }

        bool is_hubbard_pbc() const {
            return m_on_site_only_0022 && m_nnp_only_1111;
        }

        size_t nci() const {
            return ci_utils::fermion_dim(nsite(), nelec());
        }

        const size_t &nsite() const {
            return m_nsite;
        }

        bool spin_conserving_1e() const {
            return m_spin_conserving_1e;
        }

        bool spin_conserving_2e() const {
            return m_spin_conserving_2e;
        }

        bool spin_conserving() const {
            return m_spin_conserving_1e && m_spin_conserving_2e;
        }

        const size_t &nelec() const {
            return m_nelec;
        }

        const bool &complex_valued() const {
            return m_complex_valued;
        }

        const size_t &int_2e_rank() const {
            return m_int_2e_rank;
        }

        buffered::FrmOnv guess_reference(const int &spin_level) const;


        /*
        void foreach_connection(const fields::Onv<0> &src_onv, const Terms::body_fn_t &body,
                                bool get_h, bool h_nonzero_only, bool include_diagonal) const {
            m_terms->foreach_connection(src_onv, body, get_h, h_nonzero_only, include_diagonal);
        }

        void foreach_connection(const fields::Onv<1> &src_onv, const Terms::body_fn_t &body,
                                bool get_h, bool h_nonzero_only, bool include_diagonal) const {
            m_terms->foreach_connection(src_onv.m_frm, body, get_h, h_nonzero_only, include_diagonal);
        }
         */

        /**
         * set the referenced ONV object to the assumed Hartree--Fock determinant within the given spin sector
         * @param onv
         *  target onv object
         * @param spin
         *  spin (MS) number
         */
        void set_hf_onv(fields::FrmOnv& onv, int spin) const {
            auto nalpha = ci_utils::nalpha(nelec(), spin);
            auto nbeta = ci_utils::nbeta(nelec(), spin);
            DEBUG_ASSERT_EQ(nalpha+nbeta, nelec(), "inconsistent na, nb, nelec");
            onv.zero();
            for (size_t i=0ul; i<nalpha; ++i) onv.set({0, i});
            for (size_t i=0ul; i<nbeta; ++i) onv.set({1, i});
        }
    };

    struct Boson {
        const size_t m_nmode;
        std::vector<defs::ham_t> m_omegas;

        Boson(const size_t& nmode, const defs::ham_t& omega): m_nmode(nmode), m_omegas(nmode, omega){}

        defs::ham_t get_element_0(const fields::BosOnv &onv) const;

        defs::ham_comp_t get_energy(const fields::BosOnv &onv) const;

        defs::ham_t get_element_0(const conn::Boson &conn) const;

        defs::ham_t get_element_0(const conn::Antisym<1> &conn) const;

        defs::ham_t get_element(const conn::Boson &bonvconn) const;
        // TODO: read in from FCIDUMP-like file
        //Boson(std::string fname)

    };

    struct Coupling {
        const size_t m_nmode, m_nboson_cutoff;
        const defs::ham_t m_v;

        Coupling(size_t nmode, size_t nboson_cutoff, defs::ham_t v);

        defs::ham_t v(const size_t &p, const size_t &q, const size_t &n) const;

        defs::ham_t get_element_1(const size_t& p, const size_t& imode, const size_t& com) const;

        defs::ham_t get_element_1(const conn::Antisym<1> &aconn) const;

        defs::ham_t get_element(const conn::Antisym<1> &aconn) const;
    };
};


#endif
#endif //M7_HAMILTONIANPARTS_H
