//
// Created by jhalson on 06/04/2021.
//

#ifndef M7_STATICTWF_H
#define M7_STATICTWF_H


#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include <src/core/parallel/Reducible.h>
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"

struct StaticTwf {
    std::vector<defs::ham_t> m_numerator;
    std::vector<defs::ham_t> m_numerator_total;
    std::vector<defs::ham_t> m_denominator;
    std::vector<defs::ham_t> m_denominator_total;
    double_t m_fermion_double_occ_penalty_factor;
    double_t m_boson_occ_penalty_factor;
    size_t m_nsite;

    StaticTwf(size_t npart, size_t nsite, double_t fermion_factor=0.0, double_t boson_factor=0.0);

    void add(const Hamiltonian<0> &ham, const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight, const fields::Onv<0> &onv);

    void add(const Hamiltonian<1> &ham, const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight, const fields::Onv<1> &onv);

    defs::ham_t evaluate_static_twf(const fields::Onv<0> &onv) const;

    defs::ham_t evaluate_static_twf(const fields::Onv<1> &onv) const;

    void reduce();


};




#endif //M7_STATICTWF_H
