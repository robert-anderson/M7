//
// Created by rja on 01/06/2021.
//

#ifndef M7_FOREACHCONNECTION_H
#define M7_FOREACHCONNECTION_H


#include <src/core/connection/Connections.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/hamiltonian/Hamiltonian.h>

#include <utility>
#include <src/core/basis/AbelianGroup.h>
#include <src/core/basis/Suites.h>


namespace foreach_conn2 {

    typedef std::function<void(const conn::FrmOnv&)> frm_fn_t;
    typedef std::function<void(const conn::FrmOnv&, defs::ham_t)> frm_h_fn_t;
    typedef std::function<void(const conn::BosOnv&)> bos_fn_t;
    typedef std::function<void(const conn::BosOnv&, defs::ham_t)> bos_h_fn_t;
    typedef std::function<void(const conn::FrmBosOnv&)> frmbos_fn_t;
    typedef std::function<void(const conn::FrmBosOnv&, defs::ham_t)> frmbos_h_fn_t;

    struct Base {
        const Hamiltonian &m_ham;
        OccupiedOrbitals m_occ;
        VacantOrbitals m_vac;
        suite::Conns m_conns;

        Base(const Hamiltonian &ham) :
                m_ham(ham), m_occ(ham.nsite()), m_vac(ham.nsite()), m_conns(ham.nsite()){}

        virtual void operator()(const fields::FrmOnv& mbf, const frm_fn_t& body_fn) = 0;
        virtual void operator()(const fields::BosOnv& mbf, const bos_fn_t& body_fn) = 0;
        virtual void operator()(const fields::FrmBosOnv& mbf, const frmbos_fn_t& body_fn) = 0;

        /*
         * the following methods simply adapt the virtual methods for computation of the matrix elements
         */
        void operator()(const fields::FrmOnv& mbf, const frm_h_fn_t& body_fn, bool nonzero_h_only) {
            frm_fn_t fn = [&](const conn::FrmOnv& conn) {
                auto helem = m_ham.get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                body_fn(conn, helem);
            };
            (*this)(mbf, fn);
        };
        void operator()(const fields::BosOnv& mbf, const bos_h_fn_t& body_fn, bool nonzero_h_only) {
            bos_fn_t fn = [&](const conn::BosOnv& conn) {
                auto helem = m_ham.get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                body_fn(conn, helem);
            };
            (*this)(mbf, fn);
        };
        void operator()(const fields::FrmBosOnv& mbf, const frmbos_h_fn_t& body_fn, bool nonzero_h_only) {
            frmbos_fn_t fn = [&](const conn::FrmBosOnv& conn) {
                auto helem = m_ham.get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                body_fn(conn, helem);
            };
            (*this)(mbf, fn);
        };
    };

    namespace frm {
        /**
         * no symmetry is respected in this basic case
         */
        struct Fermion : Base {
            using Base::operator();

            Fermion(const Hamiltonian &ham) : Base(ham) {}

        protected:
            virtual void singles(conn::FrmOnv &conn, const std::function<void()> &fn) const {
                for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
                    for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
                        conn.clear();
                        conn.add(m_occ[iocc], m_vac[ivac]);
                        fn();
                    }
                }
            }

            virtual void singles(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) const {
                singles(conn, fn);
            }

            virtual void doubles(conn::FrmOnv &conn, const std::function<void()> &fn) const {
                for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
                    for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
                        for (size_t jocc = 0ul; jocc < iocc; ++jocc) {
                            for (size_t jvac = 0ul; jvac < ivac; ++jvac) {
                                conn.clear();
                                conn.add(m_occ[jocc], m_occ[iocc], m_vac[jvac], m_vac[ivac]);
                                fn();
                            }
                        }
                    }
                }
            }

            virtual void operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
                m_occ.update(mbf);
                m_vac.update(mbf);
                singles(conn, fn);
                doubles(conn, fn);
            }

        public:
            void operator()(const fields::FrmOnv &mbf, const frm_fn_t &body_fn) override {
                auto fn = [&]() { body_fn(m_conns.m_frmonv); };
                (*this)(mbf, m_conns.m_frmonv, fn);
            }

            void operator()(const fields::BosOnv &mbf, const bos_fn_t &body_fn) override {}

            void operator()(const fields::FrmBosOnv &mbf, const frmbos_fn_t &body_fn) override {
                auto fn = [&]() { body_fn(m_conns.m_frmbosonv); };
                (*this)(mbf.m_frm, m_conns.m_frmbosonv.m_frm, fn);
            }
        };

        /**
         * make use of spin and point-group symmetry in the loop over connections
         */
        struct SpinSym : Fermion {
            using Fermion::operator();

            SpinSym(const Hamiltonian &ham) : Fermion(ham) {}

        protected:
            void singles(conn::FrmOnv &conn, const std::function<void()> &fn) const override {
                Fermion::singles(conn, fn);
            }

            void doubles(conn::FrmOnv &conn, const std::function<void()> &fn) const override {
                Fermion::doubles(conn, fn);
            }
        };

        struct Hubbard1D : Fermion {
            using Fermion::operator();
            const bool m_pbc;

            Hubbard1D(const Hamiltonian &ham) : Fermion(ham), m_pbc(ham.m_frm.m_int_1(0, ham.nsite() - 1) != 0.0) {}

        protected:
            void
            singles(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) const override {
                for (const auto &occ: m_occ.inds()) {
                    size_t neighbor;
                    neighbor = conn_utils::left(occ, m_ham.nsite(), m_pbc);
                    if (neighbor != ~0ul && !mbf.get(neighbor)) {
                        conn.set(occ, neighbor);
                        fn();
                    }
                    neighbor = conn_utils::right(occ, m_ham.nsite(), m_pbc);
                    if (neighbor != ~0ul && !mbf.get(neighbor)) {
                        conn.set(occ, neighbor);
                        fn();
                    }
                }
            }

        public:

            void operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) override {
                m_occ.update(mbf);
                singles(mbf, conn, fn);
            }
        };
    }
}


namespace foreach_conn {

    struct Base {
        const Hamiltonian &m_ham;
        std::function<void(defs::ham_t)> m_body_fn;
        const bool m_get_h, m_nonzero_h_only;
        defs::ham_t m_helem = 0.0;

        Base(const Hamiltonian &ham, std::function<void(defs::ham_t)> body_fn, bool get_h = true, bool nonzero_h_only = true) :
                m_ham(ham), m_body_fn(std::move(body_fn)), m_get_h(get_h), m_nonzero_h_only(nonzero_h_only) {}
    };

    struct Fermion : Base {
        conn::FrmOnv &m_conn;
        const size_t m_nsite;
        OccupiedOrbitals m_occ;
        VacantOrbitals m_vac;

        Fermion(const Hamiltonian &ham, conn::FrmOnv &conn,
                std::function<void(defs::ham_t)> body_fn, bool get_h = true, bool nonzero_h_only = true) :
                Base(ham, body_fn, get_h, nonzero_h_only),
                m_conn(conn), m_nsite(ham.nsite()), m_occ(m_nsite), m_vac(m_nsite) {
        }

    protected:
        bool update_helem(const fields::FrmOnv &src_onv) {
            if (m_get_h) m_helem = m_ham.m_frm.get_element(src_onv, m_conn);
            return !m_nonzero_h_only || !consts::float_is_zero(m_helem);
        }

        void body_fn(const fields::FrmOnv &src_onv) {
            if (update_helem(src_onv)) m_body_fn(m_helem);
        }

    public:
        virtual void operator()(const fields::FrmOnv &src_onv) {
            ASSERT(!src_onv.is_zero());
            m_occ.update(src_onv);
            m_vac.update(src_onv);

            for (auto &iocc: m_occ.inds()) {
                for (auto &ivac: m_vac.inds()) {
                    // singles
                    m_conn.set(iocc, ivac);
                    body_fn(src_onv);
                    for (auto &jocc: m_occ.inds()) {
                        // doubles
                        if (jocc <= iocc) continue;
                        for (auto &jvac: m_vac.inds()) {
                            if (jvac <= ivac) continue;
                            m_conn.set(iocc, jocc, ivac, jvac);
                            body_fn(src_onv);
                        }
                    }
                }
            }
        }
    };

    struct Hubbard1D : Fermion {
        const bool m_pbc;

        Hubbard1D(const Hamiltonian &ham, conn::FrmOnv &conn,
                  std::function<void(defs::ham_t)> body_fn, bool pbc) :
                Fermion(ham, conn, body_fn, true, true), m_pbc(pbc) {}

        void operator()(const fields::FrmOnv &src_onv) override {
            ASSERT(!src_onv.is_zero());
            m_occ.update(src_onv);
            for (auto &iocc: m_occ.inds()) {
                auto neighbor = conn_utils::left(iocc, m_ham.nsite(), m_pbc);
                if (neighbor != ~0ul && !src_onv.get(neighbor)) {
                    // there is an orbital to the left, and it is unoccupied
                    m_conn.set(iocc, neighbor);
                    body_fn(src_onv);
                }
                neighbor = conn_utils::right(iocc, m_ham.nsite(), m_pbc);
                if (neighbor != ~0ul && !src_onv.get(neighbor)) {
                    // there is an orbital to the right, and it is unoccupied
                    m_conn.set(iocc, neighbor);
                    body_fn(src_onv);
                }
            }
        }
    };
}

#if 0
    /**
     * currently implements single changes in boson occupation for occupied fermion sites
     * e.g. Hubbard--Holstein
     */
    struct FermiBos : Base {
        const Hamiltonian<1> &m_ham;
        conn::Antisym<1> &m_conn;
        const size_t m_nsite;
        Fermion m_frm;

        FermiBos(const Hamiltonian<1> &ham, conn::Antisym<1> &conn, Fermion &&frm,
                 std::function<void(defs::ham_t)> body_fn, bool get_h = true, bool nonzero_h_only = true) :
                Base(std::move(body_fn), get_h, nonzero_h_only),
                m_ham(ham), m_conn(conn), m_nsite(ham.nsite()), m_frm(frm) {}

    protected:
        bool update_helem() {
            if (m_get_h) m_helem = m_ham.get_element(m_conn);
            return !m_nonzero_h_only || !consts::float_is_zero(m_helem);
        }

        void body_fn() {
            if (update_helem()) m_body_fn(m_helem);
        }

    public:
        void operator()(const fields::Onv<1> &src_onv) {
            m_conn.m_bonvconn.zero();
            m_frm(src_onv.m_frm);
            static_cast<FermionOnvConnection&>(m_conn).connect(src_onv.m_frm, src_onv.m_frm);
            // m_occ has already been updated by the loop over fermionic sites
            for (auto &iocc: m_frm.m_occ.inds()) {
                auto isite = iocc < m_nsite ? iocc : iocc-m_nsite;
                if (src_onv.m_bos[isite]>0) {
                    // we have a boson-annihilating coupling to the electron density
                    m_conn.m_bonvconn.zero();
                    m_conn.m_bonvconn.add(isite, -1);
                    m_conn.apply(src_onv);
                    body_fn();
                }
                if (src_onv.m_bos[isite]<m_ham.nboson_cutoff()) {
                    // we have a boson-creating coupling to the electron density
                    m_conn.m_bonvconn.zero();
                    m_conn.m_bonvconn.add(isite, 1);
                    m_conn.apply(src_onv);
                    body_fn();
                }
            }
        }
    };
}
#endif

#if 0
    /**
     * only generates connections which conserve spin and the abelian symmetry provided
     */
    struct FermionSymm : Fermion {
        const AbelianGroup m_grp;

        template<size_t nind>
        struct ForeachData {
            const AbelianGroup &m_grp;
            foreach::ctnd::Ordered<nind, false, true> m_foreach_cre_grp;
            foreach::ctnd::Ordered<nind, false, true> m_foreach_ann_grp;
            foreach::ctnd::Unrestricted<nind> m_foreach_cre_spin;
            foreach::ctnd::Unrestricted<nind> m_foreach_ann_spin;
            foreach::ctnd::body_fn_t m_body_fn;

            ForeachData(const AbelianGroup &grp, foreach::ctnd::cr_body_fn_t body_fn) :
                    m_grp(grp),
                    m_foreach_cre_grp(grp.nirrep()), m_foreach_ann_grp(grp.nirrep()),
                    m_foreach_cre_spin(2), m_foreach_ann_spin(2), m_body_fn(body_fn) {}

            bool is_conservative() const {
                if (!m_grp.is_conservative(m_foreach_ann_grp.inds(), m_foreach_cre_grp.inds()))
                    return false;
                return m_foreach_ann_spin.sum() == m_foreach_cre_spin.sum();
            }

            void operator()() {
                foreach::ctnd::chain(m_body_fn, m_foreach_cre_grp, m_foreach_ann_grp, m_foreach_cre_spin,
                                     m_foreach_ann_spin);
            }
        };


        struct DoublesForeachData : ForeachData<2> {
            DoublesForeachData(const AbelianGroup &grp) :
                    ForeachData<2>(grp, [&]() {
                        if (!is_conservative()) return;
                    }) {}
        };

        SinglesForeachData m_foreach_singles;
        DoublesForeachData m_foreach_doubles;

    public:
        FermionSymm(AbelianGroup grp, const Hamiltonian<0> &ham, conn::Antisym<0> &conn,
                    std::function<void(defs::ham_t)> body_fn, bool get_h = true, bool nonzero_h_only = true) :
                Fermion(ham, conn, body_fn, get_h, nonzero_h_only),
                m_grp(std::move(grp)),
                m_foreach_singles(grp), m_foreach_doubles(grp) {}

    public:
        void operator()(const fields::Onv<0> &src_onv)
        override;
    };

    struct Hubbard1D : Fermion {
        Hubbard1D(const Hamiltonian<0> &ham, conn::Antisym<0> &conn,
                  std::function<void(defs::ham_t)> body_fn, bool get_h = true, bool nonzero_h_only = true) :
                Fermion(ham, conn, body_fn, get_h, nonzero_h_only) {}

        void operator()(const fields::Onv<0> &src_onv) override {
            ASSERT(!src_onv.is_zero());
            m_occ.update(src_onv);
        }
    };


}

#endif

#endif //M7_FOREACHCONNECTION_H
