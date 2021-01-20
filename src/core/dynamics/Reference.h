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

class Reference : public elements::Onv<>, ra::Onv::Dynamic {
    Wavefunction &m_wf;
    const Hamiltonian<> &m_ham;
    size_t m_irow;
    size_t m_irank;

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
    NdArray<defs::wf_t, defs::ndim_wf> m_weight;
    SingleReducible<defs::wf_comp_t, defs::ndim_wf> m_candidate_weight;

public:
    Reference(const Options &m_opts, Wavefunction &wf, const Hamiltonian<> &ham, views::Onv<> &onv);

    void add_row(const size_t& irow);

    using elements::Onv<>::operator=;
    using elements::Onv<>::mpi_bcast;
    void change(const size_t& irow, const size_t& irank);

    void log_candidate_weight(const size_t& irow, const defs::wf_comp_t& candidate_weight);

    void reset();

    void reduce();

    const size_t &irow();

    const bool &in_redefinition_cycle();

    bool is_mine() const;

    bool is_connected(const views::Onv<> &onv) const;

    void add_to_numerator(const views::Onv<> &onv, const defs::wf_t &weight);

    ReductionMember<defs::wf_comp_t, defs::ndim_wf> &nwalker_at_doubles();

    ReductionMember<defs::wf_comp_t, defs::ndim_wf> &candidate_weight();

    defs::ham_t proj_energy_num() const;
    defs::ham_comp_t weight() const;
    defs::ham_comp_t proj_energy() const;

private:


    void on_row_send_(size_t irow) override;

    void on_row_recv_(size_t irow) override;
};

#endif //M7_REFERENCE_H
