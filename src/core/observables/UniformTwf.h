//
// Created by rja on 18/03/2021.
//

#ifndef M7_UNIFORMTWF_H
#define M7_UNIFORMTWF_H


#include <src/core/parallel/Reducible.h>
#include "src/core/wavefunction/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "src/core/observables/SignProblemFreeTwf.h"
#include "src/core/basis/Suites.h"

template<typename fn_t, bool enable_bosons = defs::enable_bosons>
struct ConnectionEnumerator {
    suite::Conns m_conn;
    suite::Mbfs m_mbf;
    OccupiedOrbitals m_occ;
    VacantOrbitals m_vac;
    const fn_t &m_fn;

    ConnectionEnumerator(size_t nsite, const fn_t &fn) :
            m_conn(nsite), m_mbf(nsite), m_occ(nsite), m_vac(nsite), m_fn(fn) {}


    void execute(const fields::FrmBosOnv &onv) {
        m_occ.update(onv.m_frm);
        m_vac.update(onv.m_frm);
        fn_t fn;
        auto& conn = m_conn[onv];
        auto& mbf = m_mbf[onv];
        for (auto &iocc: m_occ.inds()) {
            for (auto &ivac: m_vac.inds()) {
                conn.clear();
                conn.m_frm.add(iocc, ivac);
                conn.m_frm.apply(onv.m_frm, mbf.m_frm);
                fn(m_conn);
                for (auto &jocc: m_occ.inds()) {
                    for (auto &jvac: m_vac.inds()) {
                        conn.clear();
                        conn.m_frm.add(iocc, jocc, ivac, jvac);
                        conn.m_frm.apply(onv.m_frm, mbf.m_frm);
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

#if 0
    void add(const Hamiltonian<1> &ham,
             const fields::Numbers<defs::wf_t,
             defs::ndim_wf> &weight,
             const fields::Onv<1> &onv) override;
#endif

    void reduce() override;
};


#endif //M7_UNIFORMTWF_H
