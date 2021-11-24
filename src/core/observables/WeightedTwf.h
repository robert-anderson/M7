//
// Created by jhalson on 06/04/2021.
//

#ifndef M7_WEIGHTEDTWF_H
#define M7_WEIGHTEDTWF_H

#include "src/core/wavefunction/WalkerTable.h"
#include "src/core/observables/SignProblemFreeTwf.h"

class WeightedTwf : public SpfTwfBase{
protected:
    double m_frm_doub_occ_penalty_factor;
    double m_bos_occ_penalty_factor;

public:
    WeightedTwf(const Hamiltonian& ham, size_t npart, size_t nsite, double_t fermion_factor=0.0, double_t boson_factor=0.0);

    virtual ~WeightedTwf(){}

private:
    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight, defs::ham_t helem_sum, defs::ham_t diag_fac);

public:
    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::FrmOnv &onv) override;

    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::FrmBosOnv &onv) override;

    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::BosOnv &onv) override;

    defs::ham_t evaluate_static_twf(const field::FrmOnv &onv) const;

    defs::ham_t evaluate_static_twf(const field::FrmBosOnv &onv) const;

};




#endif //M7_WEIGHTEDTWF_H
