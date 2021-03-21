//
// Created by rja on 18/03/2021.
//

#ifndef M7_UNIFORMTWF_H
#define M7_UNIFORMTWF_H


#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include <src/core/parallel/Reducible.h>
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"

template<typename fn_t, bool enable_bosons = defs::enable_bosons>
struct ConnectionEnumerator {
    conn::Antisym<enable_bosons> m_conn;
    buffered::Onv<enable_bosons> m_work_onv;
    OccupiedOrbitals m_occ;
    VacantOrbitals m_vac;
    const fn_t &m_fn;

    ConnectionEnumerator(size_t nsite, const fn_t &fn) :
            m_conn(nsite), m_work_onv(nsite), m_occ(nsite), m_vac(nsite), m_fn(fn) {}


    void execute(const fields::Onv<1> &onv) {
        m_occ.update(onv.m_fonv);
        m_vac.update(onv.m_fonv);
        fn_t fn;
        for (auto &iocc: m_occ.inds()) {
            for (auto &ivac: m_vac.inds()) {
                m_conn.zero();
                m_conn.add(iocc, ivac);
                m_conn.apply(onv, m_work_onv);
                fn(m_conn);
                for (auto &jocc: m_occ.inds()) {
                    for (auto &jvac: m_vac.inds()) {
                        m_conn.zero();
                        m_conn.add(iocc, jocc, ivac, jvac);
                        m_conn.apply(onv, m_work_onv);
                        fn(m_conn);
                    }
                }
            }
        }
    }

};

struct UniformTwf {
    std::vector<defs::ham_t> m_numerator;
    std::vector<defs::ham_t> m_numerator_total;
    size_t m_nsite;

    UniformTwf(size_t npart, size_t nsite) :
            m_numerator(npart, 0.0), m_numerator_total(npart, 0.0), m_nsite(nsite) {}

    void add(const Hamiltonian<0> &ham, const fields::Vector<defs::wf_t> &weight, const fields::Onv<0> &onv) {
        conn::Antisym<0> conn(m_nsite);
        buffered::Onv<0> work_onv(m_nsite);
        OccupiedOrbitals occ(m_nsite);
        VacantOrbitals vac(m_nsite);

        defs::ham_t helem_sum = 0.0;
        occ.update(onv);
        vac.update(onv);

        conn.connect(onv, onv);
        helem_sum += ham.get_element(conn);
        for (auto &iocc: occ.inds()) {
            for (auto &ivac: vac.inds()) {
                conn.zero();
                conn.add(iocc, ivac);
                conn.apply(onv, work_onv);
                helem_sum += ham.get_element(conn);
                for (auto &jocc: occ.inds()) {
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
            m_numerator[ipart] += weight(ipart) * helem_sum;
        }
    }

    void add(const Hamiltonian<1> &ham, const fields::Vector<defs::wf_t> &weight, const fields::Onv<1> &onv) {
        conn::Antisym<1> conn(m_nsite);
        buffered::Onv<1> work_onv(m_nsite);
        OccupiedOrbitals occ(m_nsite);
        VacantOrbitals vac(m_nsite);

        defs::ham_t helem_sum = 0.0;
        occ.update(onv);
        vac.update(onv);

        conn.connect(onv, onv);
        helem_sum += std::abs(ham.get_element(conn));
        for (auto &iocc: occ.inds()) {

            for (int change = -1; change <= 1; change += 2) {
                const size_t imode = iocc < ham.nsite() ? iocc : iocc - ham.nsite();
                auto vacd_minus = (onv.m_bonv(imode) == 0) && (change < 0);
                auto occd_plus = (onv.m_bonv(imode) == ham.nboson_cutoff()) && (change > 0);
                if (!vacd_minus && !occd_plus){
                    conn.zero();
                    conn.m_bonvconn.add(imode, change);
                    work_onv.zero();
                    conn.apply(onv, work_onv);
                    auto com = work_onv.m_bonv(imode);
                    if (change < 0) com += change;
                    helem_sum += std::abs(ham.bc().get_element_1(imode, imode, com));
                }
            }

            for (auto &ivac: vac.inds()) {
                conn.zero();
                conn.add(iocc, ivac);
                conn.apply(onv, work_onv);
                helem_sum += std::abs(ham.get_element(conn));
//                for (auto &jocc: occ.inds()) {
//                    if (jocc <= iocc) continue;
//                    for (auto &jvac: vac.inds()) {
//                        if (jvac<=ivac) continue;
//                        conn.zero();
//                        conn.add(iocc, jocc, ivac, jvac);
//                        ASSERT(conn.ncre()==2 && conn.nann()==2)
//                        work_onv = onv;
//                        conn.apply(onv, work_onv);
//                        helem_sum += std::abs(ham.get_element(conn));
//                    }
//                }
            }
        }
        for (size_t ipart = 0ul; ipart < m_numerator.size(); ++ipart) {
            m_numerator[ipart] -= std::abs(weight(ipart)) * helem_sum;
        }
    }

    void reduce() {
        mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
        m_numerator.assign(m_numerator.size(), 0.0);
    }


};


#endif //M7_UNIFORMTWF_H
