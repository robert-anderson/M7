//
// Created by Robert J. Anderson on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H

#include <M7_lib/bilinear/Bilinears.h>
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
     * The hole/particle subspace (Dpert):
     * D is the deterministic (walker) subspace, Dpert is the union of (all perturbers acting on all members of D)
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
     * we'll call this object Dpert_maps
     *
     * Dpert_maps = empty map
     * Dpert = empty set
     * for all perturbers q:
     *   add a new element to dpert_maps with key q and value of an empty list of integer pairs
     *   for all d in D:
     *     if q acting on d destroys the state:
     *       loop
     *     if qd is not in Dpert:
     *       // we have found a new basis state
     *       insert qd into Dpert
     *     // else we are revisiting a previously generated element of the union
     *     append (index of d in D, index of qd in Dpert) pair to Dpert_maps[q]
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
     * now all setup has been set out, all that remains is to specify the actual computation of the inner product:
     *
     * let Ci be the gathered walker weights of both replicas across D
     * let Vwork be a working vector with the same length as Dpert
     * for all perturbers q:
     *   // remember that just like in RDMs, the bra and ket coefficients must be uncorrelated
     *   for ireplica in [0, 1]:
     *     ireplica_other = !ireplica
     *     Vwork = 0
     *     for all elements (ici, ipert) in Dpert_maps[q]:
     *       Vwork[ipert] += Ci[ici][ireplica]
     *     // now the (N+/-1)-electron ket has been prepared
     *     for ih in (0, 1, ..., n): // where n is the order of the moment
     *       // perform sparse, distributed multiplication - including MPI allgather
     *       Vwork = Hpert * Vwork
     *     // Vwork now contains an (N+/-1)-electron state which has been propagated n times
     *     // perform final contractions
     *     for all perturbers p:
     *       inner_product = 0
     *       for all elements (ici, ipert) in Dpert_maps[p]:
     *         inner_product += Vwork[ipert] * Ci[ici][ireplica_other]
     *     make contribution inner_product to An(+/-)[p, q]
     */
    struct FrmOpPerturbed {
        /**
         * true if the perturber is of the "hole" kind, as opposed to "particle"
         */
        const bool m_hole;
        /**
         * mapped list of basis (N+1)- or (N-1)- electron functions to aid in the setup of sparse H
         */
        typedef SingleFieldRow<field::Mbf> basis_row_t;
        typedef buffered::MappedTable<basis_row_t> basis_table_t;
        basis_table_t m_basis_table;
        /**
         * all the non-zero Hamiltonian connections and their matrix elements for the locally-stored rows (with respect to
         * columns distributed across all ranks)
         */
        sparse::dynamic::Matrix<ham_t> m_ham_matrix;
        /**
         * there is no available optimised distribution for the (N+1)- or (N-1)- electron basis so create an arbitrary
         * one so that parallel multiplication can be employed
         */
        Distribution m_dist;
        /**
         * spin-orbital indices required for p and q
         */
        const uintv_t m_peturber_inds;
        FrmOpPerturbed(sys::Basis basis, bool hole, uintv_t peturber_inds, uint_t iroot):
            m_hole(hole),
            m_basis_table(logging::format("Perturbed ({}) deterministic subspace for root {}",
              hole ? "hole" : "particle", iroot), {basis_row_t(basis, "mbf")}),
            m_dist(mpi::nrank()*6),
            m_peturber_inds(std::move(peturber_inds)){}

        /**
         * @return
         *  true if the given spin orbital has the correct occupation/vacancy status with respect to the perturber kind
         */
        bool has_required_occ(const field::FrmOnv& onv, uint_t ispinorb) {
            return m_hole ? onv.get(ispinorb) : !onv.get(ispinorb);
        }
        bool has_required_occ(const field::BosOnv&, uint_t ispinorb) {
            return false;
        }
        bool has_required_occ(const field::FrmBosOnv& onv, uint_t ispinorb) {
            return has_required_occ(onv.m_frm, ispinorb);
        }

        /**
         * if hole kind, clear the ispinorb bit, else set it
         */
        void modify_occ(field::FrmOnv& onv, uint_t ispinorb) {
            DEBUG_ASSERT_TRUE(has_required_occ(onv), ispinorb);
            if (m_hole) onv.clr(ispinorb);
            else onv.set(ispinorb);
        }
        void modify_occ(field::BosOnv& onv, uint_t ispinorb) {}
        bool modify_occ(field::FrmBosOnv& onv, uint_t ispinorb) {
            modify_occ(onv.m_frm, ispinorb);
        }

        void make_connections(const MappedTable<Walker>& walker_space, const Hamiltonian& ham, uint_t ispinorb) {
            buffered::Mbf work_mbf(walker_space.m_row.m_mbf);
            auto walker = walker_space.m_row;
            for (walker.restart(); walker; ++walker) {
                // ignore this N-electron state if it would be destroyed by application of the perturber
                if (!has_required_occ(walker.m_mbf, ispinorb)) continue;
                work_mbf = walker.m_mbf;
                modify_occ(work_mbf, ispinorb);
                // see if this perturber has already been generated
                auto& lookup = m_basis_table.lookup(work_mbf);
                if (lookup) {

                }
                else {
                    m_basis_table.insert(work_mbf);
                }
            }
        }

        void make_connections(const MappedTable<Walker>& walker_space, const Hamiltonian& ham) {
            for (uint_t ispinorb=0ul; ispinorb<m_basis.m_frm.m_nspinorb; ++ispinorb){
                make_connections(walker_space, ham, ispinorb);
            }
        }
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
        Wavefunction& m_wf;
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

        Subspace(const conf::Semistochastic& opts, Wavefunction& wf, uint_t iroot);

        virtual ~Subspace() {}

        /**
         * add a row to the subspace on this MPI rank (not collectively), and reflect this status in the flags of the WF row
         * @param row
         *  row of the wavefunction store which is pointing at the MBF to add into the subspace
         */
        void add_(Walker& row);

        void select_highest_weighted();

        void select_l1_norm_fraction();

        void make_connections(const Hamiltonian& ham, const Bilinears& bilinears);

        void make_rdm_contribs(Rdms& rdms, const shared_rows::Walker* hf);

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
        void init(const Hamiltonian& ham, const Bilinears& bilinears, Wavefunction& wf, uint_t icycle);

        void update();

        void project(double tau);

        void make_rdm_contribs(Rdms& rdms, const shared_rows::Walker* hf);
    };
}

#endif //M7_DETERMINISTICSUBSPACE_H
