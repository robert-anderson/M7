//
// Created by rja on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include "src/core/parallel/RankAllocator.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/nd/NdArray.h"
#include "WalkerTable.h"
#include "Wavefunction.h"

class Reference : public WalkerMappedTable::DynamicRow {
    const Hamiltonian<> &m_ham;
    const Wavefunction &m_wf;

    mutable conn::Antisym<> m_aconn;

    /*
     * If a candidate for redefinition of the reference is found, then
     * its weight and row within m_list must be stored
     */
    size_t m_irow_candidate;
    const double m_redefinition_thresh;
    bool m_redefinition_cycle;

    ReductionSyndicate m_summables;
    ReductionMember<defs::ham_t, defs::ndim_wf> m_proj_energy_num;
    ReductionMember<defs::wf_comp_t, defs::ndim_wf> m_nwalker_at_doubles;
    SingleReducible<defs::wf_comp_t, defs::ndim_wf> m_candidate_weight;

public:
    Reference(const Options &m_opts, const Hamiltonian<> &ham,
              const Wavefunction &wf, Table::Loc loc);

    const views::Onv<> get_onv() const {
        return m_ac.m_onv(0);
    }

    const defs::wf_t& get_weight(const size_t& iroot, const size_t& ireplica) const {
        return m_ac.m_weight(0, iroot, ireplica);
    }

    void add_row(const size_t& irow);

    //void change(const size_t& irow, const size_t& irank);

    //void log_candidate_weight(const size_t& irow, const defs::wf_comp_t& candidate_weight);

    void reset();

    void reduce();

    //const bool &in_redefinition_cycle();

    bool is_connected(const views::Onv<> &onv) const;

    void add_to_numerator(const views::Onv<> &onv, const defs::wf_t &weight);

    ReductionMember<defs::wf_comp_t, defs::ndim_wf> &nwalker_at_doubles();

    ReductionMember<defs::wf_comp_t, defs::ndim_wf> &candidate_weight();

    defs::ham_t proj_energy_num() const;
    defs::ham_comp_t weight() const;
    defs::ham_comp_t proj_energy() const;

};

#endif //M7_REFERENCE_H
