//
// Created by rja on 18/03/2021.
//

#include "UniformTwf.h"

UniformTwf::UniformTwf(const Hamiltonian& ham, size_t npart, size_t nsite) :
        SpfTwfBase(ham, npart, nsite){}

void UniformTwf::add(const Numbers<defs::wf_t, defs::ndim_wf> &weight, defs::ham_t helem_sum) {
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += weight[ipart] * helem_sum;
    }
}

void UniformTwf::add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::FrmOnv &onv) {
    defs::ham_t helem_sum = m_ham.get_element(onv);
    auto fn = [&helem_sum](defs::ham_t helem){helem_sum+=helem;};
    m_foreach_conn->foreach(onv, fn);
    add(weight, helem_sum);
}

void UniformTwf::add(const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::FrmBosOnv &onv) {
    defs::ham_t helem_sum = m_ham.get_element(onv);
    auto fn = [&helem_sum](defs::ham_t helem){helem_sum+=helem;};
    m_foreach_conn->foreach(onv, fn);
    add(weight, helem_sum);
}

void UniformTwf::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
}