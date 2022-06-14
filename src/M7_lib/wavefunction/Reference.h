//
// Created by Robert J. Anderson on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include <M7_lib/observables/RefExcits.h>
#include <M7_lib/basis/Suites.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/parallel/ReductionMember.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/nd/NdArray.h>

#include "WalkerTable.h"
#include "Wavefunction.h"

class Reference : public Wavefunction::SharedRow {
    const Hamiltonian &m_ham;
    const Wavefunction &m_wf;
    /**
     * index to the "part" of the wavefunction for which this object tracks the reference row
     */
    const size_t m_ipart;

    /**
     * work space for computing connections to other ONVs via the Hamiltonian
     */
    mutable suite::Conns m_conn;

    /**
     * If a candidate for redefinition of the reference is found, then its weight and row within m_list must be stored
     */
    size_t m_irow_candidate;
    /**
     * candidate weights are averaged over replicas so that replicas on the same root will redefine to the same row
     */
    defs::wf_comp_t m_candidate_abs_weight = 0.0;
    /**
     * default scale factor defining when the candidate weight has grown to the magnitude at which it must be made the
     * new reference for this WF part
     */
    const double m_redefinition_thresh;

    ReductionSyndicate m_summables;
    Reduction<defs::ham_t> m_proj_energy_num;
    Reduction<defs::wf_comp_t> m_nwalker_at_doubles;

public:
    Reference(const conf::Reference &opts, const Hamiltonian &ham,
              const Wavefunction &wf, size_t ipart, TableBase::Loc loc);

    const field::Mbf& get_mbf() const;

    size_t occupied_ncycle(const size_t& icycle) const {
        return m_global.m_row.occupied_ncycle(icycle);
    }

    /**
     * Perform a special loop over stored rows of the WF, updating the reference connection flags of each non-zero row.
     */
    void update_ref_conn_flags();

    /**
     * If the candidate with the highest absolute weight across all MPI ranks has a weight larger than that of the current
     * reference by a specified factor, then accept the candidate and let both the dynamic row set and reference
     * connection flags in the WF m_store table reflect this change
     * @param redefinition_thresh
     *  scale factor by which the candidate weight must exceed that of the current ref
     */
    void accept_candidate(double redefinition_thresh = 0.0);

    /**
     * add contributions from the current m_wf.m_store.m_row
     */
    void contrib_row();

    /**
     * reset variables to begin a fresh MC cycle
     */
    void begin_cycle();

    void end_cycle();

    /**
     * @param mbf
     *  reference to MBF
     * @return
     *  true if the matrix element of the Hamiltonian between the current reference and the argument is non-zero
     */
    bool is_connected(const field::Mbf &mbf) const;
    /**
     * @param mbf
     *  reference to MBF
     * @return
     *  excitation signature from ref to arg
     */
    size_t exsig(const field::Mbf &mbf) const;

    /**
     * occupied ONVs connected to the reference must contribute to the numerator inner product <ref | H | onv>
     * @param onv
     * @param weights
     */
    void make_numerator_contribs(const field::Mbf &onv, const defs::wf_t& weight);

    const defs::wf_comp_t &nwalker_at_doubles();

    const defs::ham_t& proj_energy_num() const;

    const defs::wf_t& weight() const;

    /**
     * this method includes the current weight in the average, bringing the normalized average up to date.
     * @param icycle
     *  cycle on which average is being used
     * @param ipart
     *  wf part index
     * @return
     *  normalized average weight
     */
    defs::wf_t norm_average_weight(const size_t& icycle, const size_t& ipart) const;

    const WalkerTableRow& row() const {
        return m_global.m_row;
    }
};

/**
 * Wavefunctions have many parts in general, and each of these may require different reference ONVs, so here we define
 * a vector of the above-defined single part Reference class
 */
struct References {
    std::vector<Reference> m_refs;
    buffered::Numbers<defs::ham_t, defs::ndim_wf> m_proj_energy_nums;
    buffered::Numbers<defs::wf_t, defs::ndim_wf> m_weights;

    References(const conf::Reference &opts, const Hamiltonian &ham, const Wavefunction &wf, std::vector<TableBase::Loc> locs);

    const Reference& operator[](const size_t& ipart) const;

    void begin_cycle();

    void end_cycle();

    void contrib_row();

    std::vector<bool> is_connected(const field::Mbf &onv) const;

    const field::Numbers<defs::ham_t, defs::ndim_wf>& proj_energy_nums();

    const field::Numbers<defs::wf_t, defs::ndim_wf>& weights();

};

#endif //M7_REFERENCE_H
