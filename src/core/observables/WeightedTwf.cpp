//
// Created by jhalson on 06/04/2021.
//

#include "WeightedTwf.h"

WeightedTwf::WeightedTwf(size_t npart, size_t nsite, double_t fermion_factor, double_t boson_factor) :
        SpfTwfBase(npart, nsite),
        m_fermion_double_occ_penalty_factor(fermion_factor),
        m_boson_occ_penalty_factor(boson_factor)
        {}

void WeightedTwf::add(const Hamiltonian<0> &ham, const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::FrmOnv &onv) {
    conn::FrmOnv conn(m_nsite);
    buffered::FrmOnv work_onv(m_nsite);
    OccupiedOrbitals occ(m_nsite);
    VacantOrbitals vac(m_nsite);

    defs::ham_t helem_sum = 0.0;
    occ.update(onv);
    vac.update(onv);

    auto this_twf = evaluate_static_twf(onv);
    helem_sum += this_twf*ham.get_element(onv);
    for (auto &iocc: occ.inds()) {
        for (auto &ivac: vac.inds()) {
            conn.clear();
            conn.add(iocc, ivac);
            helem_sum += evaluate_static_twf(work_onv)*ham.get_element(onv, conn);
            for (auto &jocc: occ.inds()) {
                if (jocc <= iocc) continue;
                for (auto &jvac: vac.inds()) {
                    if (jvac<=ivac) continue;
                    conn.clear();
                    conn.add(iocc, jocc, ivac, jvac);
                    helem_sum += evaluate_static_twf(work_onv)*ham.get_element(onv, conn);
                }
            }
        }
    }
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += weight[ipart] * helem_sum;
        m_denominator[ipart] += weight[ipart] * this_twf;
    }
}

defs::ham_t WeightedTwf::evaluate_static_twf(const fields::Onv<0> &onv) const{
    size_t num_double_occ_sites = 0;
    for(size_t isite=0; isite < m_nsite; isite++) {
        num_double_occ_sites += (onv.get({0, isite}) and onv.get({1, isite}));
    }
    return std::exp(-m_fermion_double_occ_penalty_factor*num_double_occ_sites);
}

defs::ham_t WeightedTwf::evaluate_static_twf(const fields::Onv<1> &onv) const{
    size_t total_boson_occ = 0;
    for(size_t isite=0; isite < m_nsite; isite++) {
        total_boson_occ += onv.m_bos[isite];
    }
    return std::exp(-m_boson_occ_penalty_factor*total_boson_occ) + evaluate_static_twf(onv.m_frm);
}

#if 0
void WeightedTwf::add(const Hamiltonian<1> &ham, const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::Onv<1> &onv) {
    conn::Antisym<1> conn(m_nsite);
    buffered::Onv<1> work_onv(m_nsite);
    OccupiedOrbitals occ(m_nsite);
    VacantOrbitals vac(m_nsite);

    defs::ham_t helem_sum = 0.0;
    defs::ham_t helem = 0.0;
    occ.update(onv);
    vac.update(onv);

    conn.connect(onv, onv);

    helem = ham.get_element_0(conn);
    auto this_twf = evaluate_static_twf(onv);
    helem_sum += this_twf * helem;


    for (auto &iocc: occ.inds()) {
        const size_t imode = iocc < ham.nsite() ? iocc : iocc - ham.nsite();
        for (int change = -1; change <= 1; change += 2) {
            auto vacd_minus = (onv.m_bos[imode] == 0) && (change < 0);
            auto occd_plus = (onv.m_bos[imode] == ham.nboson_cutoff()) && (change > 0);
            if (!vacd_minus && !occd_plus){
                auto com = onv.m_bos[imode];
                if (change<0) com+=change;

                work_onv.zero();
                conn.zero();
                conn.m_bonvconn.add(imode, change);
                conn.apply(onv, work_onv);

                helem = ham.bc().get_element_1(imode, imode, com);
                helem_sum-=std::abs(helem)*evaluate_static_twf(work_onv);
            }
        }

        for (auto &ivac: vac.inds()) {
            conn.zero();
            conn.add(iocc, ivac);
            work_onv.zero();
            conn.apply(onv, work_onv);
            helem = ham.get_element_1(conn);
            helem_sum -= std::abs(helem)*evaluate_static_twf(work_onv);
        }
    }
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += std::abs(weight[ipart]) * helem_sum;
        m_denominator[ipart] += this_twf * std::abs(weight[ipart]);
    }
}
#endif