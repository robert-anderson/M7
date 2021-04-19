//
// Created by jhalson on 06/04/2021.
//

#ifndef M7_STATICTWF_H
#define M7_STATICTWF_H


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

struct StaticTwf {
    std::vector<defs::ham_t> m_numerator;
    std::vector<defs::ham_t> m_numerator_total;
    std::vector<defs::ham_t> m_denominator;
    std::vector<defs::ham_t> m_denominator_total;
    double_t m_fermion_double_occ_penalty_factor;
    double_t m_boson_occ_penalty_factor;
    size_t m_nsite;

    StaticTwf(size_t npart, size_t nsite, double_t fermion_factor=0.0, double_t boson_factor=0.0);

    void add(const Hamiltonian<0> &ham, const fields::Vector<defs::wf_t> &weight, const fields::Onv<0> &onv);

    void add(const Hamiltonian<1> &ham, const fields::Vector<defs::wf_t> &weight, const fields::Onv<1> &onv);

    defs::ham_t evaluate_static_twf(const fields::Onv<0> &onv) const;

    defs::ham_t evaluate_static_twf(const fields::Onv<1> &onv) const;

    void reduce();


};




#endif //M7_STATICTWF_H
