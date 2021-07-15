//
// Created by rja on 18/03/2021.
//

#ifndef M7_UNIFORMTWF_H
#define M7_UNIFORMTWF_H


#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include <src/core/parallel/Reducible.h>
#include "src/core/wavefunction/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/observables/SignProblemFreeTwf.h"

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
        m_occ.update(onv.m_frm);
        m_vac.update(onv.m_frm);
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

class UniformTwf : public SpfTwfBase {
public:
    UniformTwf(size_t npart, size_t nsite);

    virtual ~UniformTwf(){}

    void add(const Hamiltonian<0> &ham,
             const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
             const fields::Onv<0> &onv) override;

    void add(const Hamiltonian<1> &ham,
             const fields::Numbers<defs::wf_t,
             defs::ndim_wf> &weight,
             const fields::Onv<1> &onv) override;

    void reduce() override;
};


#endif //M7_UNIFORMTWF_H
