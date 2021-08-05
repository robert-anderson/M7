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

    void add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const fields::FrmOnv &onv) override;

    void add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const fields::FrmBosOnv &onv) override;


    defs::ham_t evaluate_static_twf(const fields::FrmOnv &onv) const;

    defs::ham_t evaluate_static_twf(const fields::FrmBosOnv &onv) const;

};




#endif //M7_WEIGHTEDTWF_H
