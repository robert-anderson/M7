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

#if 0
/**
 * captures only the data from Walker relevant to semistochastic propagation
 */
struct DeterministicDataRow : Row {
    field::Mbf m_mbf;
    field::Numbers<wf_t, c_ndim_wf> m_weight;

    field::Mbf &key_field();

    DeterministicDataRow(const Wavefunction &wf);

    static void load_fn(const Walker &source, DeterministicDataRow &local);
};
#endif

/**
 * the implementation of the semistochastic adaptation calls for a subspace of many-body basis functions to be
 * designated as a deterministic subspace within which the exact propagator is applied to update the CI coefficients
 */
struct DeterministicSubspace : shared_rows::Set<Walker> {
    /**
     * options from the configuration object
     */
    const conf::Semistochastic &m_opts;
    /**
     * the deterministic subspace is just a selection of rows from an MPI-distributed wavefunction object
     */
    typedef shared_rows::Set<Walker> base_t;
    /**
     * need to maintain a reference to the wavefunction so the affected coefficients can be gathered and updated
     */
    Wavefunction &m_wf;
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

    void make_rdm_contrib(Rdms &rdms, const shared_rows::Walker *hf, const sparse::Element& elem);


public:

    DeterministicSubspace(const conf::Semistochastic &opts, Wavefunction &wf, uint_t iroot);

    virtual ~DeterministicSubspace() {}

    /**
     * add a row to the subspace on this MPI rank (not collectively), and reflect this status in the flags of the WF row
     * @param row
     *  row of the wavefunction store which is pointing at the MBF to add into the subspace
     */
    void add_(Walker &row);

    void select_highest_weighted();

    void select_l1_norm_fraction();

    void make_connections(const Hamiltonian &ham, const Bilinears &bilinears);

    void make_rdm_contribs(Rdms &rdms, const shared_rows::Walker* hf);

    /**
      * for every deterministically-propagated row on this MPI rank, update its value.
      * @param tau
      *  timestep
      */
    void project(double tau);

    void save(const hdf5::NodeWriter& nw) const {
        m_all.m_row.m_mbf.save(nw,  logging::format("root_{}", m_iroot), true);
    }

    void load(const hdf5::NodeReader& /*nr*/) const {
//        using namespace hdf5::dataset;
//        buffered::Table<SingleFieldRow<Mbf>> load_table("detsub loader", SingleFieldRow<Mbf>(m_all.m_row.m_mbf));
//        /*
//         *
//        typedef std::function<buf_t*(const ListFormat& format, uint_t max_nitem_per_op)> load_prep_fn;
//        typedef std::function<void(const buf_t* src, uint_t nitem)> load_fill_fn;
//         */
//        auto load_prep_fn = [&](const ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
//
//        };
//        auto load_fill_fn = [&](const ListFormat& format, uint_t max_nitem_per_op) -> buf_t* {
//
//        };
//
//
//        nr.load_dataset(logging::format("root_{}", m_iroot), hdf5::dataset::load_prep_fn);
//        load_table.m_row.m_field.load()
//        load_table.load()
//        m_all.m_row.m_mbf.save(nw,  logging::format("root_{}", m_iroot), true);
    }
};

struct DeterministicSubspaces {
    /**
     * options from the configuration object
     */
    const conf::Semistochastic &m_opts;
    /**
     * epoch for all deterministic subspaces. not used beyond the provision of logging
     */
    Epoch m_epoch;
    /**
     * one instance of DeterministicSubspace per root
     */
    v_t<std::unique_ptr<DeterministicSubspace>> m_detsubs;

    DeterministicSubspaces(const conf::Semistochastic &opts);

    operator bool() const;

    /**
     * create subspaces, perform the relevant selections, and call make_connections on each
     */
    void init(const Hamiltonian &ham, const Bilinears &bilinears, Wavefunction &wf, uint_t icycle);

    void update();

    void project(double tau);

    void make_rdm_contribs(Rdms &rdms, const shared_rows::Walker *hf);
};

#endif //M7_DETERMINISTICSUBSPACE_H
