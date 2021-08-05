//
// Created by rja on 18/03/2021.
//

#ifndef M7_UNIFORMTWF_H
#define M7_UNIFORMTWF_H


#include "src/core/wavefunction/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/observables/SignProblemFreeTwf.h"
#include "src/core/basis/Suites.h"

class UniformTwf : public SpfTwfBase {
public:
    UniformTwf(const Hamiltonian &ham, size_t npart, size_t nsite);

    virtual ~UniformTwf() {}

    void add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const fields::FrmOnv &onv) override;

    void add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const fields::FrmBosOnv &onv) override;

    void reduce() override;
};


#endif //M7_UNIFORMTWF_H
