//
// Created by rja on 18/03/2021.
//

#include "UniformTwf.h"

UniformTwf::UniformTwf(const Hamiltonian& ham, size_t npart, size_t nsite) :
        SpfTwfBase(ham, npart, nsite){}

void UniformTwf::add(const Numbers<defs::wf_t, defs::ndim_wf> &weight, defs::ham_t helem_sum) {
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += std::abs(weight[ipart]) * helem_sum;
    }
}

void UniformTwf::add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const field::FrmOnv &onv) {
    defs::ham_t helem_sum = m_ham.get_element(onv);
    auto fn = [&helem_sum](defs::ham_t helem){helem_sum-=std::abs(helem);};
    m_excit_iters.foreach<field::FrmOnv>(onv, fn);
    add(weight, helem_sum);
}

void UniformTwf::add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const field::FrmBosOnv &onv) {
    defs::ham_t helem_sum = m_ham.get_element(onv);
    auto fn = [&helem_sum](const conn::FrmBosOnv& conn, defs::ham_t helem){
//        std::cout << conn.m_frm.ann() << " -> " << conn.m_frm.cre() << "    "
//        << conn.m_bos.m_ann.to_vector() << " -> " << conn.m_bos.m_cre.to_vector() << std::endl;
        helem_sum-=std::abs(helem);
    };
//    std::cout << onv << std::endl;
    m_excit_iters.foreach<field::FrmBosOnv>(onv, fn, false);
    add(weight, helem_sum);
}

void UniformTwf::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
}