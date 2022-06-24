//
// Created by Robert J. Anderson on 18/03/2021.
//

#include "UniformTwf.h"

UniformTwf::UniformTwf(const Hamiltonian& ham, uint_t npart) : SpfTwfBase(ham, npart){}

void UniformTwf::add(const Numbers<wf_t, c_ndim_wf> &weight, ham_t helem_sum) {
    for (uint_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += std::abs(weight[ipart]) * helem_sum;
    }
}

void UniformTwf::add(const field::Numbers<wf_t, c_ndim_wf> &weight,
                     const field::FrmOnv &onv) {
    ham_t helem_sum = m_ham.get_element(onv);
    conn::FrmOnv conn(m_ham.m_basis.size());
    auto fn = [&](){
        auto helem = m_ham.get_element(onv, conn);
        helem_sum-=std::abs(helem);
    };
    m_conn_iters.loop(conn, onv, fn);
    add(weight, helem_sum);
}

void UniformTwf::add(const field::Numbers<wf_t, c_ndim_wf> &weight,
                     const field::FrmBosOnv &onv) {
    ham_t helem_sum = m_ham.get_element(onv);
    conn::FrmBosOnv conn(m_ham.m_basis.size());
    auto fn = [&](){
        auto helem = m_ham.get_element(onv, conn);
        helem_sum-=std::abs(helem);
    };
    m_conn_iters.loop(conn, onv, fn);
    add(weight, helem_sum);
}

void UniformTwf::add(const field::Numbers<wf_t, c_ndim_wf> &/*weight*/,
                     const field::BosOnv &/*onv*/) {
    ABORT("not yet implemented");
}

void UniformTwf::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
}