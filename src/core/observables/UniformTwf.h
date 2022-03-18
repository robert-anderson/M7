//
// Created by rja on 18/03/2021.
//

#ifndef M7_UNIFORMTWF_H
#define M7_UNIFORMTWF_H


#include "wavefunction/WalkerTable.h"
#include "hamiltonian/Hamiltonian.h"
#include "observables/SignProblemFreeTwf.h"
#include "basis/Suites.h"

class UniformTwf : public SpfTwfBase {
public:
    UniformTwf(const Hamiltonian &ham, size_t npart);

    virtual ~UniformTwf() {}

private:
    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight, defs::ham_t helem_sum);

public:
    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::FrmOnv &onv) override;

    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::FrmBosOnv &onv) override;

    void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const field::BosOnv &onv) override;

    void reduce() override;
};


#endif //M7_UNIFORMTWF_H
