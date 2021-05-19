//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/nd/NdArray.h"
#include "WalkerTable.h"
#include "Wavefunction.h"

class Reference : public Wavefunction::SharedRow {
    const Hamiltonian<> &m_ham;
    const Wavefunction &m_wf;
    /*
     * index to the "part" of the wavefunction for which this object tracks the reference row
     */
    const size_t m_ipart;

    mutable conn::Antisym<> m_aconn;

    /*
     * If a candidate for redefinition of the reference is found, then
     * its weight and row within m_list must be stored
     */
    size_t m_irow_candidate;
    defs::wf_t m_candidate_abs_weight = 0.0;
    const double m_redefinition_thresh;
    bool m_redefinition_cycle;

    ReductionSyndicate m_summables;
    NdReduction<defs::ham_t, defs::ndim_wf> m_proj_energy_num;
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nwalker_at_doubles;

public:
    Reference(const Options &m_opts, const Hamiltonian<> &ham,
              const Wavefunction &wf, size_t ipart, TableBase::Loc loc);

    const fields::Onv<>& get_onv() const {
        return m_global.m_row.m_onv;
    }

    void accept_candidate(double redefinition_thresh = 0.0) {
        std::vector<defs::wf_t> gather(mpi::nrank());
        mpi::all_gather(m_candidate_abs_weight, gather);
        MPI_ASSERT(m_candidate_abs_weight==gather[mpi::irank()], "Gather error");
        size_t irank = std::distance(gather.begin(), std::max_element(gather.begin(), gather.end()));
        mpi::bcast(m_irow_candidate, irank);
        if (gather[irank] > std::abs(weight()[irank]*redefinition_thresh)){
            log::debug("Changing the reference ONV. current weight: {}, candidate: {}", weight()[irank], gather[irank]);
            change({irank, m_irow_candidate});
            //MPI_ASSERT(std::abs(m_wf.m_store.m_weight(m_irow_candidate, 0, 0))==m_candidate_abs_weight, "");
        }
        ASSERT(std::abs(weight()[irank])==m_candidate_abs_weight);
        m_candidate_abs_weight = 0.0;
    }

    void update() override;

    /**
     * add contributions from the current row of Wavefunction::m_store
     */
    void add_row();

    //void change(const size_t& irow, const size_t& irank);

    //void log_candidate_weight(const size_t& irow, const defs::wf_comp_t& candidate_weight);

    void begin_cycle();

    void end_cycle();

    //const bool &in_redefinition_cycle();

    bool is_connected(const fields::Onv<> &onv) const;

    void add_to_numerator(const fields::Onv<> &onv, const fields::Numbers<defs::ham_t, defs::ndim_wf> &weights);

    NdReduction<defs::wf_comp_t, defs::ndim_wf> &nwalker_at_doubles();

    NdReduction<defs::wf_comp_t, defs::ndim_wf> &candidate_weight();

    const fields::Numbers<defs::ham_t, defs::ndim_wf>& proj_energy_num() const;

    const fields::Numbers<defs::ham_t, defs::ndim_wf>& weight() const;

    const fields::Numbers<defs::ham_t, defs::ndim_wf>& average_weight() const;

    const WalkerTableRow& row() const {
        return m_global.m_row;
    }

};

#endif //M7_REFERENCE_H
