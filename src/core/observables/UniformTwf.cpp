//
// Created by rja on 18/03/2021.
//

#include "UniformTwf.h"

UniformTwf::UniformTwf(size_t npart, size_t nsite) :
        m_numerator(npart, 0.0), m_numerator_total(npart, 0.0),
        m_nsite(nsite){}

void UniformTwf::add(const Hamiltonian<0> &ham,
                     const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::Onv<0> &onv) {
    conn::Antisym<0> conn(m_nsite);
    buffered::Onv<0> work_onv(m_nsite);
    OccupiedOrbitals occ(m_nsite);
    VacantOrbitals vac(m_nsite);

    defs::ham_t helem_sum = 0.0;
    occ.update(onv);
    vac.update(onv);

    // diagonal
    conn.connect(onv, onv);
    helem_sum += ham.get_element(conn);
    for (auto &iocc: occ.inds()) {
        for (auto &ivac: vac.inds()) {
            // singles
            conn.zero();
            conn.add(iocc, ivac);
            conn.apply(onv, work_onv);
            helem_sum += ham.get_element(conn);
            for (auto &jocc: occ.inds()) {
                // doubles
                if (jocc <= iocc) continue;
                for (auto &jvac: vac.inds()) {
                    if (jvac<=ivac) continue;
                    conn.zero();
                    conn.add(iocc, jocc, ivac, jvac);
                    conn.apply(onv, work_onv);
                    helem_sum += ham.get_element(conn);
                }
            }
        }
    }
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += weight[ipart] * helem_sum;
    }
}

void UniformTwf::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
}


void UniformTwf::add(const Hamiltonian<1> &ham, const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
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
    helem_sum += helem;

    for (auto &iocc: occ.inds()) {
        const size_t imode = iocc < ham.nsite() ? iocc : iocc - ham.nsite();
        for (int change = -1; change <= 1; change += 2) {
            auto vacd_minus = (onv.m_bos[imode] == 0) && (change < 0);
            auto occd_plus = (onv.m_bos[imode] == ham.nboson_cutoff()) && (change > 0);
            if (!vacd_minus && !occd_plus){
                auto com = onv.m_bos[imode];
                if (change<0) com+=change;
                helem = ham.bc().get_element_1(imode, imode, com);
                helem_sum-=std::abs(helem);
            }
        }

        for (auto &ivac: vac.inds()) {
            conn.zero();
            conn.add(iocc, ivac);
            work_onv.zero();
            conn.apply(onv, work_onv);
            helem = ham.get_element_1(conn);
            helem_sum -= std::abs(helem);
        }
    }
    for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
        m_numerator[ipart] += std::abs(weight[ipart]) * helem_sum;
    }
}
