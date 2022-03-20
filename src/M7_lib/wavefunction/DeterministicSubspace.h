//
// Created by rja on 16/06/2020.
//

#ifndef M7_DETERMINISTICSUBSPACE_H
#define M7_DETERMINISTICSUBSPACE_H

#include <M7_lib/bilinear/Bilinears.h>
#include <M7_lib/linalg/Sparse.h>
#include <M7_lib/hamiltonian/FrmHam.h>
#include <M7_lib/field/Fields.h>

#include "WalkerTable.h"
#include "Wavefunction.h"

/**
 * captures only the data from WalkerTableRow relevant to semistochastic propagation
 */
struct DeterministicDataRow : Row {
    field::Mbf m_mbf;
    field::Numbers<defs::wf_t, defs::ndim_wf> m_weight;

    field::Mbf &key_field();

    DeterministicDataRow(const Wavefunction &wf);

    static void load_fn(const WalkerTableRow &source, DeterministicDataRow &local);
};

/**
 * the implementation of the semistochastic adaptation calls for a subspace of many-body basis functions to be
 * designated as a deterministic subspace within which the exact propagator is applied to update the CI coefficients
 */
struct DeterministicSubspace : Wavefunction::PartSharedRowSet<DeterministicDataRow> {
    /**
     * options from the configuration object
     */
    const fciqmc_config::Semistochastic &m_opts;
    /**
     * the deterministic subspace is just a selection of rows from an MPI-distributed wavefunction object
     */
    typedef Wavefunction::PartSharedRowSet<DeterministicDataRow> base_t;
    /**
     * need to maintain a reference to the wavefunction so the affected coefficients can be gathered and updated
     */
    Wavefunction &m_wf;
    /**
     * all the non-zero Hamiltonian connections and their matrix elements for the locally-stored rows (with respect to
     * columns distributed across all ranks)
     */
    sparse::Matrix<defs::ham_t> m_ham_matrix;
    /**
     * in general RDMs take contributions from MBF connections which correspond to H matrix elements of zero. These
     * are stored here. We do NOT duplicate elements from the sparse_ham here, only those MBF pairs which may
     * contribute to RDM elements, but which are H-unconnected
     */
    sparse::Network m_rdm_network;
    /**
     * associated WF root index
     */
    const size_t m_iroot;
private:
    /**
     * the part indices associated with this subspace, will contain two entries if WF uses replication, else one
     */
    const defs::inds m_iparts;

    defs::inds make_iparts();

    void make_rdm_contribs(Rdms &rdms, const field::Mbf &ref, const defs::inds& icol_list);


public:

    DeterministicSubspace(const fciqmc_config::Semistochastic &opts, Wavefunction &wf, size_t iroot);

    virtual ~DeterministicSubspace() {}

    /**
     * add a row to the subspace on this MPI rank (not collectively), and reflect this status in the flags of the WF row
     * @param row
     *  row of the wavefunction store which is pointing at the MBF to add into the subspace
     */
    void add_(WalkerTableRow &row);

    void build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears);

    void build_connections(const Hamiltonian &ham, const Bilinears &bilinears);

    void make_rdm_contribs(Rdms &rdms, const field::Mbf &ref);

    /**
      * for every deterministically-propagated row on this MPI rank, update its value.
      * @param tau
      *  timestep
      */
    void project(double tau);
};

struct DeterministicSubspaces {
    /**
     * options from the configuration object
     */
    const fciqmc_config::Semistochastic &m_opts;
    /**
     * epoch for all deterministic subspaces. not used beyond the provision of logging
     */
    Epoch m_epoch;
    /**
     * one instance of DeterministicSubspace per root
     */
    std::vector<std::unique_ptr<DeterministicSubspace>> m_detsubs;

    DeterministicSubspaces(const fciqmc_config::Semistochastic &opts);

    operator bool() const;

    void build_from_most_occupied(const Hamiltonian &ham, const Bilinears &bilinears, Wavefunction &wf, size_t icycle);

    void update();

    void project(double tau);

    void make_rdm_contribs(Rdms &rdms, const field::Mbf &ref);
};

#endif //M7_DETERMINISTICSUBSPACE_H
