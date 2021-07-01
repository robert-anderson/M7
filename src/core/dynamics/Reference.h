//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include <src/core/observables/RefExcits.h>
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/nd/NdArray.h"
#include "WalkerTable.h"
#include "Wavefunction.h"

class Reference : public Wavefunction::SharedRow {
    const Hamiltonian<> &m_ham;
    const Wavefunction &m_wf;
    /**
     * index to the "part" of the wavefunction for which this object tracks the reference row
     */
    const size_t m_ipart;

    /**
     * work space for computing connections to other ONVs via the Hamiltonian
     */
    mutable conn::Antisym<> m_aconn;

    /**
     * If a candidate for redefinition of the reference is found, then its weight and row within m_list must be stored
     */
    size_t m_irow_candidate;
    defs::wf_t m_candidate_abs_weight = 0.0;
    /**
     * default scale factor defining when the candidate weight has grown to the magnitude at which it must be made the
     * new reference for this WF part
     */
    const double m_redefinition_thresh;

    ReductionSyndicate m_summables;
    Reduction<defs::ham_t> m_proj_energy_num;
    Reduction<defs::wf_comp_t> m_nwalker_at_doubles;

public:
    Reference(const fciqmc_config::Reference &opts, const Hamiltonian<> &ham,
              const Wavefunction &wf, size_t ipart, TableBase::Loc loc);

    const fields::Onv<>& get_onv() const;

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

    /**
     * redefine to candidate ref if needed, and perform necessary MPI reductions
     */
    void end_cycle();

    /**
     * @param onv
     *  reference to ONV
     * @return
     *  true if the matrix element of the Hamiltonian between the current reference and the argument is non-zero
     */
    bool is_connected(const fields::Onv<> &onv) const;

    bool connection_phase(const fields::Onv<> &onv) const;

    /**
     * occupied ONVs connected to the reference must contribute to the numerator inner product <ref | H | onv>
     * @param onv
     * @param weights
     */
    void make_numerator_contribs(const fields::Onv<> &onv, const defs::wf_t& weight);

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

    References(const fciqmc_config::Reference &opts, const Hamiltonian<> &ham, const Wavefunction &wf, std::vector<TableBase::Loc> locs):
            m_proj_energy_nums(wf.m_format.m_shape), m_weights(wf.m_format.m_shape){
        ASSERT(locs.size()==wf.m_format.m_nelement);
        m_refs.reserve(wf.m_format.m_nelement);
        for (size_t ipart=0ul; ipart<wf.m_format.m_nelement; ++ipart) m_refs.emplace_back(opts, ham, wf, ipart, locs[ipart]);
        ASSERT(m_refs.size()==wf.npart());
    }

    const Reference& operator[](const size_t& ipart) const {
        return m_refs[ipart];
    }

    void begin_cycle() {
        for (auto& ref: m_refs) ref.begin_cycle();
    }

    void end_cycle() {
        for (auto& ref: m_refs) ref.end_cycle();
    }

    void contrib_row() {
        for (auto& ref: m_refs) ref.contrib_row();
    }

    std::vector<bool> is_connected(const fields::Onv<>& onv) const {
        std::vector<bool> out;
        out.reserve(m_refs.size());
        for (size_t ipart=0ul; ipart<m_refs.size(); ++ipart)
            out.push_back(m_refs[ipart].is_connected(onv));
        return out;
    }

    const fields::Numbers<defs::ham_t, defs::ndim_wf>& proj_energy_nums() {
        size_t ipart = 0ul;
        for (auto& ref: m_refs) m_proj_energy_nums[ipart++] = ref.proj_energy_num();
        return m_proj_energy_nums;
    }

    const fields::Numbers<defs::wf_t, defs::ndim_wf>& weights() {
        size_t ipart = 0ul;
        for (auto& ref: m_refs) m_weights[ipart++] = ref.weight();
        return m_weights;
    }

};

#endif //M7_REFERENCE_H
