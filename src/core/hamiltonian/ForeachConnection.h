//
// Created by rja on 01/06/2021.
//

#ifndef M7_FOREACHCONNECTION_H
#define M7_FOREACHCONNECTION_H


#include <src/core/connection/Connections.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <utility>
#include <src/core/basis/AbelianGroup.h>
#include <src/core/basis/Suites.h>

struct Hamiltonian;

namespace foreach_conn {

    typedef std::function<void(const conn::FrmOnv&)> frm_fn_t;
    typedef std::function<void(const conn::FrmOnv&, defs::ham_t)> frm_h_fn_t;
    typedef std::function<void(const conn::FrmOnv&, const fields::FrmOnv&)> frm_d_fn_t;
    typedef std::function<void(const conn::FrmOnv&, const fields::FrmOnv&, defs::ham_t)> frm_dh_fn_t;

    typedef std::function<void(const conn::BosOnv&)> bos_fn_t;
    typedef std::function<void(const conn::BosOnv&, defs::ham_t)> bos_h_fn_t;
    typedef std::function<void(const conn::BosOnv&, const fields::BosOnv&)> bos_d_fn_t;
    typedef std::function<void(const conn::BosOnv&, const fields::BosOnv&, defs::ham_t)> bos_dh_fn_t;

    typedef std::function<void(const conn::FrmBosOnv&)> frmbos_fn_t;
    typedef std::function<void(const conn::FrmBosOnv&, defs::ham_t)> frmbos_h_fn_t;
    typedef std::function<void(const conn::FrmBosOnv&, const fields::FrmBosOnv&)> frmbos_d_fn_t;
    typedef std::function<void(const conn::FrmBosOnv&, const fields::FrmBosOnv&, defs::ham_t)> frmbos_dh_fn_t;

    struct Base {
        const Hamiltonian &m_ham;
        OccupiedOrbitals m_occ;
        VacantOrbitals m_vac;
        suite::Conns m_conns;
        suite::Mbfs m_mbfs;

        explicit Base(const Hamiltonian &ham);

        virtual void operator()(const fields::FrmOnv& mbf, const frm_fn_t& body_fn) = 0;
        virtual void operator()(const fields::BosOnv& mbf, const bos_fn_t& body_fn) = 0;
        virtual void operator()(const fields::FrmBosOnv& mbf, const frmbos_fn_t& body_fn) = 0;

        /*
         * adapt the above virtual methods for computation of the matrix elements
         */
        void operator()(const fields::FrmOnv& mbf, const frm_h_fn_t& body_fn, bool nonzero_h_only);
        void operator()(const fields::BosOnv& mbf, const bos_h_fn_t& body_fn, bool nonzero_h_only);
        void operator()(const fields::FrmBosOnv& mbf, const frmbos_h_fn_t& body_fn, bool nonzero_h_only);
        /*
         * adapt the above virtual methods for computation of the connected MBF
         */
        void operator()(const fields::FrmOnv& mbf, const frm_d_fn_t& body_fn) {
            auto& dst_mbf = m_mbfs[mbf];
            frm_fn_t fn = [&](const conn::FrmOnv &conn) {
                conn.apply(mbf, dst_mbf);
                body_fn(conn, dst_mbf);
            };
            (*this)(mbf, fn);
        }
        void operator()(const fields::BosOnv& mbf, const bos_d_fn_t& body_fn) {
            auto& dst_mbf = m_mbfs[mbf];
            bos_fn_t fn = [&](const conn::BosOnv &conn) {
                conn.apply(mbf, dst_mbf);
                body_fn(conn, dst_mbf);
            };
            (*this)(mbf, fn);
        }
    };

    namespace frm {
        /**
         * no symmetry is respected in this basic case
         */
        struct Fermion : Base {
            using Base::operator();

            explicit Fermion(const Hamiltonian &ham) : Base(ham) {}

        protected:
            virtual void singles(conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void singles(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void doubles(conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn);

        public:
            void operator()(const fields::FrmOnv &mbf, const frm_fn_t &body_fn) override;

            void operator()(const fields::BosOnv &mbf, const bos_fn_t &body_fn) override {}

            void operator()(const fields::FrmBosOnv &mbf, const frmbos_fn_t &body_fn) override;
        };

        /**
         * make use of spin and point-group symmetry in the loop over connections
         */
        struct SpinSym : Fermion {
            using Fermion::operator();

            explicit SpinSym(const Hamiltonian &ham) : Fermion(ham) {}

        protected:
            void singles(conn::FrmOnv &conn, const std::function<void()> &fn) override;

            void doubles(conn::FrmOnv &conn, const std::function<void()> &fn) override;
        };

        struct Hubbard1D : Fermion {
            using Fermion::operator();
            const bool m_pbc;

            explicit Hubbard1D(const Hamiltonian &ham);

        protected:
            void singles(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) override;

        public:

            void operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) override;
        };
    }
}

#endif //M7_FOREACHCONNECTION_H
