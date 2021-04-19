//
// Created by jhalson on 06/04/2021.
//

#include "StaticTwf.h"

StaticTwf::StaticTwf(size_t npart, size_t nsite, double_t fermion_factor, double_t boson_factor) :
        m_numerator(npart, 0.0), m_numerator_total(npart, 0.0),
        m_denominator(npart, 0.0), m_denominator_total(npart, 0.0),
        m_fermion_double_occ_penalty_factor(fermion_factor),
        m_boson_occ_penalty_factor(boson_factor),
        m_nsite(nsite){}

void StaticTwf::add(const Hamiltonian<0> &ham, const fields::Vector<defs::wf_t> &weight,
                     const fields::Onv<0> &onv) {
    conn::Antisym<0> conn(m_nsite);
    buffered::Onv<0> work_onv(m_nsite);
    OccupiedOrbitals occ(m_nsite);
    VacantOrbitals vac(m_nsite);

    defs::ham_t helem_sum = 0.0;
    occ.update(onv);
    vac.update(onv);

    conn.connect(onv, onv);
    auto this_twf = evaluate_static_twf(onv);
    helem_sum += this_twf*ham.get_element(conn);
    for (auto &iocc: occ.inds()) {
        for (auto &ivac: vac.inds()) {
            conn.zero();
            conn.add(iocc, ivac);
            conn.apply(onv, work_onv);
            helem_sum += evaluate_static_twf(work_onv)*ham.get_element(conn);
            for (auto &jocc: occ.inds()) {
                if (jocc <= iocc) continue;
                for (auto &jvac: vac.inds()) {
                    if (jvac<=ivac) continue;
                    conn.zero();
                    conn.add(iocc, jocc, ivac, jvac);
                    conn.apply(onv, work_onv);
                    helem_sum += evaluate_static_twf(work_onv)*ham.get_element(conn);
                }
            }
        }
    }
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += weight(ipart) * helem_sum;
        m_denominator[ipart] += weight(ipart) * this_twf;
    }
}

void StaticTwf::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
}

defs::ham_t StaticTwf::evaluate_static_twf(const fields::Onv<0> &onv) const{
    size_t num_double_occ_sites = 0;
    for(size_t isite=0; isite < m_nsite; isite++) {
        num_double_occ_sites += (onv.get(0, isite) and onv.get(1, isite));
    }
    return std::exp(-m_fermion_double_occ_penalty_factor*num_double_occ_sites);
}

defs::ham_t StaticTwf::evaluate_static_twf(const fields::Onv<1> &onv) const{
    size_t total_boson_occ = 0;
    for(size_t isite=0; isite < m_nsite; isite++) {
        total_boson_occ += onv.m_bonv(isite);
    }
    return std::exp(-m_boson_occ_penalty_factor*total_boson_occ) + evaluate_static_twf(onv.m_fonv);
}

void StaticTwf::add(const Hamiltonian<1> &ham, const fields::Vector<defs::wf_t> &weight,
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
            auto vacd_minus = (onv.m_bonv(imode) == 0) && (change < 0);
            auto occd_plus = (onv.m_bonv(imode) == ham.nboson_cutoff()) && (change > 0);
            if (!vacd_minus && !occd_plus){
                auto com = onv.m_bonv(imode);
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
        m_numerator[ipart] += std::abs(weight(ipart)) * helem_sum;
        m_denominator[ipart] += this_twf * std::abs(weight(ipart))
    }
}
