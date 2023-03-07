//
// Created by Robert J. Anderson on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H

#include "M7_lib/mae/Maes.h"
#include "M7_lib/linalg/sparse/Sparse.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"
#include <M7_lib/field/Fields.h>

#include "WalkerTable.h"
#include "Wavefunction.h"

namespace deterministic {

    /**
     * Spectral moments call for the preparation of the perturbed deterministic subspace. this class initializes
     * objects required for the deterministic non-diagonal contributions
     *
     * parallel sparse matrix x vector multiplications of the form Mv = u are done as follows:
     * if v and u are column-vectors with n and m elements respectively, then M is an (m x n) matrix
     * v must be available in full on all processes, so its inner product with rows of M can be evaluated
     * elements of u, and rows of M can be distributed
     *
     * for the walker subspace, we only need to gather the instantaneous value of walkers in the subspace from all ranks
     * then update the locally stored portion of the walker table.
     *
     * here, we must generate a non-propagating N+1 or N-1 electron subspace, in which no walkers are stored. since
     * we have two uncontracted indices (p, q) the task of computing the contributions could be cast as a sparse,
     * distributed matrix multiplication, but looping over indices p, q and performing matrix-vector multiplications
     * is preferred for the reason of hermiticity: we only need to compute the inner product for the spin-conserving
     * pairs (p, q) with p >= q.
     *
     * The hole/particle subspace (pert_basis):
     * D is the deterministic (walker) subspace, pert_basis is the union of (all perturbers acting on all members of D)
     *
     * for a given perturber index, there is a one-to-many-or-none relationship between the N-electron MBFs and the
     * (N+/-1)-electrons MBFs. This relationship must be sparsely cached in a structure e.g.:
     * {
     *   spinorb_i : [(i_ci_0, i_pert_0), (i_ci_1, i_pert_1), (i_ci_2, i_pert_2)],
     *   spinorb_j : [(j_ci_0, j_pert_0), (j_ci_1, j_pert_1), (j_ci_2, j_pert_2), (j_ci_3, j_pert_3)]
     * }
     * in this example, spinorb_i corresponds to 3 non-state-destroying applications of the perturber on D the i_ci_*
     * elements of this row can occur only once, but the i_pert_* elements can occur many times or not at all
     *
     * we'll call this object basis_maps
     *
     * basis_maps = empty map
     * pert_basis = empty set
     * for all perturbers q:
     *   add a new element to dpert_maps with key q and value of an empty list of integer pairs
     *   for all d in D:
     *     if q acting on d destroys the state:
     *       loop
     *     if qd is not in pert_basis:
     *       // we have found a new basis state
     *       insert qd into pert_basis
     *     // else we are revisiting a previously generated element of the union
     *     append (index of d in D, index of qd in pert_basis, fermi phase of perturber application) to basis_maps[q]
     *
     * this algorithm is implemented in the setup_basis method
     *
     * once this code has run, we have the perturber basis and the mapping between walker basis and perturber basis
     * initialized, now we must build Hpert, Hamiltonian projected into this perturber basis. to do this, distribute the
     * rows of Hpert evenly between processors, and build in the same way that the sparse H is built in the Subspace
     * class, but with an important difference:
     * due to the multiple applications of H, the idea of a diagonal element becomes more complicated. The normal
     * semi-stochastic space does not perform walker death, rather this is left to the propagator. However, it is
     * better to accumulate the true diagonals (all factors are of the form Hii) due to D space MBF pairs in this class
     * than to attempt to exclude them. Therefore, at the end of a block averaging cycle (or walker death, but D space
     * walkers are protected and never allowed to die) we only perform contributions to the true diagonals of the
     * moments if the walker in question is not a member of D.
     *
     * the setup of the hamiltonian in the perturbed basis is implemented in the setup_ham method
     *
     * now all setup has been set out, all that remains is to specify the actual computation of the inner product:
     *
     * let Ci be the gathered walker weights of both replicas across D
     * let Vwork be a working vector with the same length as pert_basis
     * for all perturbers q:
     *   // remember that just like in RDMs, the bra and ket coefficients must be uncorrelated
     *   for ireplica in [0, 1]:
     *     ireplica_other = !ireplica
     *     Vwork = 0
     *     for all elements (ici, ipert) in basis_maps[q]:
     *       Vwork[ipert] += Ci[ici][ireplica]
     *     // now the (N+/-1)-electron ket has been prepared
     *     for ih in [0, 1, ..., n]: // where n is the order of the moment
     *       // Vwork now contains an (N+/-1)-electron state which has been propagated ih times
     *       // perform final contractions
     *       for all perturbers p:
     *          inner_product = 0
     *          for all elements (ici, ipert) in basis_maps[p]:
     *            inner_product += Vwork[ipert] * Ci[ici][ireplica_other]
     *          make contribution inner_product to An(+/-)[p, q]
     *       // perform sparse, distributed multiplication - including MPI allgather
     *       Vwork = Hpert * Vwork
     *
     * this algorithm is implemented in make_contribs
     */
    struct FrmOpPerturbed {
        const sys::frm::Basis m_frm_basis;
        /**
         * true if the Hamiltonian conserves Ms2
         */
        const bool m_ms2_conserve;
        /**
         * true if the perturber is of the "hole" kind, as opposed to "particle"
         */
        const bool m_hole;
        /**
         * mapped list of basis (N+/-1)-electron functions which comprise the perturbed basis
         */
        typedef SingleFieldRow<field::Mbf> pert_basis_row_t;
        typedef buffered::MappedTable<pert_basis_row_t> pert_basis_table_t;
        pert_basis_table_t m_pert_basis_table;
        /**
         * all the non-zero Hamiltonian connections and their matrix elements for the locally-stored rows (with respect to
         * columns distributed across all ranks)
         */
        sparse::dynamic::Matrix<ham_t> m_ham_pert;
        /**
         * local offset of the first row of the distributed m_ham_pert on this MPI rank
         */
        uint_t m_pert_basis_displ;
        /**
         * spin-orbital indices required for p and q
         */
        const uintv_t m_select_pert_inds;
        /**
         * pairs of indices linking the walker basis to the pertuber basis. keys are given by m_select_pert_inds see
         * class commentary above for details
         */
        struct BasisMapEntry {
            uint_t m_iwalker;
            uint_t m_iperturbed;
            bool m_phase;
        };
        std::map<uint_t, std::list<BasisMapEntry>> m_basis_maps;
        /**
         * working vector that extends over entire perturbed space
         */
        v_t<wf_t> m_full_work_vec;
        /**
         * working vector that extends over just the part of the perturbed space for which this rank is responsible
         */
        v_t<wf_t> m_part_work_vec;

        FrmOpPerturbed(sys::Sector sector, bool hole, uintv_t select_pert_inds, uint_t iroot);

    private:
        uint_t full_basis_size() const {
            return m_pert_basis_table.nrow_in_use();
        }

        uint_t part_basis_size() const {
            return m_ham_pert.nrow();
        }

        /**
         * @param ispinorb
         *  spin orbital index of the pertuber function
         * @return
         *  true if the perturber belongs to the selected set
         */
        bool is_selected_pertuber(uint_t ispinorb) const;

        static void set_conn(conn::FrmOnv& conn, uint_t ispinorb, bool hole) {
            conn.clear();
            auto& ops = hole ? conn.m_ann : conn.m_cre;
            ops.add(ispinorb);
        }
        static void set_conn(conn::BosOnv&, uint_t, bool) {}
        static void set_conn(conn::FrmBosOnv& conn, uint_t ispinorb, bool hole) {
            set_conn(conn.m_frm, ispinorb, hole);
        }

        /**
         * initialize the perturbed basis and the maps from walker D-space to perturbed space indices
         */
        void setup_basis(const MappedTable<Walker>& walker_subspace, uint_t ispinorb);

        void setup_basis(const MappedTable<Walker>& walker_subspace) {
            for (auto& ispinorb: m_select_pert_inds) setup_basis(walker_subspace, ispinorb);
        }

        /**
         * initialize sparse Hamiltonian representing the projection of full H into the perturbed deterministic space
         */
        void setup_ham(const Hamiltonian& h);

        /**
         * rezero, then populate the perturbed space working vector
         */
        void fill_full_vec(const MappedTable<Walker>& walker_subspace, uint_t ispinorb_right, uint_t ipart_right);

        /**
         * m_full_work_vec contains the perturbed D-space vector multiplied by some power of H, contract it with the
         * D-space vector projected onto the left-perturber to obtain a moment element estimate
         */
        wf_t contract(const MappedTable<Walker>& walker_subspace, uint_t ispinorb_left, uint_t ipart_left);

        /**
         * multiply m_full_work_vec by the sparse, perturbed space Hamiltonian
         */
        void project_ham();

    public:

        void setup(const Hamiltonian& h, const MappedTable<Walker>& walker_subspace) {
            setup_basis(walker_subspace);
            setup_ham(h);
        }

        /**
         * make all contributions (diagonal and off-diagonal) to all orders of hole and particle spectral moments
         * given a pair of replica indices
         */
        void make_contribs(SpecMoms& spec_moms, const MappedTable<Walker>& walker_subspace, uint_t ipart, uint_t ipart_replica);
    };

    /**
     * the implementation of the semistochastic adaptation calls for a subspace of many-body basis functions to be
     * designated as a deterministic subspace within which the exact propagator is applied to update the CI coefficients
     */
    struct Subspace : shared_rows::Set<Walker> {
        /**
         * options from the configuration object
         */
        const conf::Semistochastic& m_opts;
        /**
         * the deterministic subspace is just a selection of rows from an MPI-distributed wavefunction object
         */
        typedef shared_rows::Set<Walker> base_t;
        /**
         * need to maintain a reference to the wavefunction so the affected coefficients can be gathered and updated
         */
        wf::Vectors& m_wf;
        /**
         * all the non-zero Hamiltonian connections and their matrix elements for the locally-stored rows (with respect to
         * columns distributed across all ranks)
         */
        sparse::dynamic::Matrix<ham_t> m_ham_matrix;
        /**
         * in general RDMs take contributions from MBF connections which correspond to H matrix elements of zero. These
         * are stored here. We do NOT duplicate elements from the sparse_ham here, only those MBF pairs which may
         * contribute to RDM elements, but which are H-unconnected
         */
        sparse::dynamic::Network m_rdm_network;

        std::unique_ptr<FrmOpPerturbed> m_frm_hole_perturbed;
        std::unique_ptr<FrmOpPerturbed> m_frm_particle_perturbed;
        /**
         * associated WF root index
         */
        const uint_t m_iroot;

        Walker m_local_row;

    private:
        /**
         * the part indices associated with this subspace, will contain two entries if WF uses replication, else one
         */
        const uintv_t m_iparts;

        uintv_t make_iparts();

        void make_rdm_contrib(Rdms& rdms, const shared_rows::Walker* hf, const sparse::Element& elem);


    public:

        Subspace(const conf::Semistochastic& opts, wf::Vectors& wf, uint_t iroot);

        virtual ~Subspace() {}

        /**
         * add a row to the subspace on this MPI rank (not collectively), and reflect this status in the flags of the WF row
         * @param row
         *  row of the wavefunction store which is pointing at the MBF to add into the subspace
         */
        void add_(Walker& row);

        void select_highest_weighted();

        void select_l1_norm_fraction();

        void make_connections(const Hamiltonian& ham, const Rdms& rdms);
        void make_connections(const Hamiltonian& ham, const SpecMoms& rdms);

        void make_rdm_contribs(Rdms& rdms, const shared_rows::Walker* hf);

        void make_spec_mom_contribs(SpecMoms& spec_moms) {
            if (!spec_moms || !spec_moms.m_accum_epoch) return;
            REQUIRE_TRUE(m_frm_hole_perturbed.get(), "uninitialized");
            REQUIRE_TRUE(m_frm_particle_perturbed.get(), "uninitialized");
            m_frm_hole_perturbed->make_contribs(spec_moms, gathered(), 0, 1);
            m_frm_hole_perturbed->make_contribs(spec_moms, gathered(), 1, 0);
            m_frm_particle_perturbed->make_contribs(spec_moms, gathered(), 0, 1);
            m_frm_particle_perturbed->make_contribs(spec_moms, gathered(), 1, 0);
        }

        /**
          * for every deterministically-propagated row on this MPI rank, update its value.
          * @param tau
          *  timestep
          */
        void project(double tau);

        void save(const hdf5::NodeWriter& nw) const {
            gathered().m_row.m_mbf.save(nw, logging::format("root_{}", m_iroot), true);
        }

    };

    struct Subspaces {
        /**
         * options from the configuration object
         */
        const conf::Semistochastic& m_opts;
        /**
         * epoch for all deterministic subspaces. not used beyond the provision of logging
         */
        Epoch m_epoch;
        /**
         * one instance of Subspace per root
         */
        v_t<std::unique_ptr<Subspace>> m_detsubs;

        Subspaces(const conf::Semistochastic& opts);

        operator bool() const;

        /**
         * create subspaces, perform the relevant selections, and call make_connections on each
         */
        void init(const Hamiltonian& ham, const Maes& maes, wf::Vectors& wf, uint_t icycle);

        void update();

        void project(double tau);

        void make_rdm_contribs(Rdms& rdms, const shared_rows::Walker* hf);

        void make_spec_mom_contribs(SpecMoms& spec_moms);
    };
}

#endif //M7_DETERMINISTICSUBSPACE_H
