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
#include "Hamiltonian.h"

using namespace fields;

namespace foreach_conn {

    template<typename mbf_t>
    using fn_c_t = std::function<void(const conn::from_field_t<mbf_t>&)>;
//    template<typename mbf_t>
//    using fn_cd_t = std::function<void(const conn::from_field_t<mbf_t>&, const mbf_t&)>;
    template<typename mbf_t>
    using fn_d_t = std::function<void(const mbf_t&)>;
    template<typename mbf_t>
    using fn_ch_t = std::function<void(const conn::from_field_t<mbf_t>&, defs::ham_t)>;
    //using fn_h_t = std::function<void(defs::ham_t)>;
    template<typename mbf_t>
    using fn_cdh_t = std::function<void(const conn::from_field_t<mbf_t>&, const mbf_t&, defs::ham_t)>;
    template<typename mbf_t>
    using fn_dh_t = std::function<void(const mbf_t&, defs::ham_t)>;

    struct Base {
        const Hamiltonian &m_ham;
        OccupiedOrbitals m_occ;
        VacantOrbitals m_vac;
        suite::Conns m_conns;
        suite::Mbfs m_mbfs;

        explicit Base(const Hamiltonian &ham);
        virtual ~Base(){}

        virtual void foreach(const FrmOnv& mbf, const fn_c_t<FrmOnv>& body_fn) = 0;
        virtual void foreach(const FrmBosOnv& mbf, const fn_c_t<FrmBosOnv>& body_fn) = 0;
        virtual void foreach(const BosOnv& mbf, const fn_c_t<BosOnv>& body_fn) = 0;

        defs::ham_t get_element(const FrmOnv& mbf, const conn::FrmOnv& conn);
        defs::ham_t get_element(const FrmBosOnv& mbf, const conn::FrmBosOnv& conn);
        defs::ham_t get_element(const BosOnv& mbf, const conn::BosOnv& conn);

        /*
         * adapt the above virtual methods for computation of the connected MBF
         */
//        template<typename mbf_t>
//        void foreach(const mbf_t& mbf, const fn_cd_t<mbf_t>& body_fn) {
//            auto& dst_mbf = m_mbfs[mbf];
//            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
//                conn.apply(mbf, dst_mbf);
//                body_fn(conn, dst_mbf);
//            };
//            this->foreach(mbf, fn);
//        }

        /*
         * adapt the above virtual methods for computation of the matrix elements
         */
//        template<typename mbf_t>
//        void foreach(const mbf_t& mbf, const fn_ch_t<mbf_t>& body_fn, bool nonzero_h_only){
//            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
//                auto helem = get_element(mbf, conn);
//                if (nonzero_h_only && consts::float_is_zero(helem)) return;
//                body_fn(conn, helem);
//            };
//            this->foreach(mbf, fn);
//        }

        /*
         * adapt the above virtual methods for computation of the matrix elements, and loop over calls to a function
         * which ONLY requires the matrix elements
         */
//        template<typename mbf_t>
//        void foreach(const mbf_t& mbf, const fn_h_t& body_fn){
//            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
//                auto helem = get_element(mbf, conn);
//                if (consts::float_is_zero(helem)) return;
//                body_fn(helem);
//            };
//            this->foreach(mbf, fn);
//        }

        /*
         * adapt the above virtual methods for computation of both the connected MBF and the matrix elements
         */
        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_cdh_t<mbf_t>& body_fn, bool nonzero_h_only){
            auto& dst_mbf = m_mbfs[mbf];
            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
                auto helem = get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                conn.apply(mbf, dst_mbf);
                body_fn(conn, dst_mbf, helem);
            };
            this->foreach(mbf, fn);
        }

        /*
         * adapt the above virtual methods for computation of both the connected MBF and the matrix elements and loop
         * over calls to a function which ONLY requires the connected MBF
         */
//        template<typename mbf_t>
//        void foreach(const mbf_t& mbf, const fn_d_t<mbf_t>& body_fn, bool nonzero_h_only){
//            auto& dst_mbf = m_mbfs[mbf];
//            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
//                auto helem = get_element(mbf, conn);
//                if (nonzero_h_only && consts::float_is_zero(helem)) return;
//                conn.apply(mbf, dst_mbf);
//                body_fn(dst_mbf);
//            };
//            this->foreach(mbf, fn);
//        }

        /*
         * adapt the above virtual methods for computation of both the connected MBF and the matrix elements and loop
         * over calls to a function which ONLY requires the connected MBF and the matrix element, but not the connection
         * itself
         */
//        template<typename mbf_t>
//        void foreach(const mbf_t& mbf, const fn_dh_t<mbf_t>& body_fn, bool nonzero_h_only){
//            auto& dst_mbf = m_mbfs[mbf];
//            auto fn = [&](const conn::from_field_t<mbf_t> &conn) {
//                auto helem = get_element(mbf, conn);
//                if (nonzero_h_only && consts::float_is_zero(helem)) return;
//                conn.apply(mbf, dst_mbf);
//                body_fn(dst_mbf, helem);
//            };
//            this->foreach(mbf, fn);
//        }
    };

    namespace frm {
        /**
         * no symmetry is respected in this basic case
         */
        struct Fermion : Base {
            using Base::foreach;

            explicit Fermion(const Hamiltonian &ham) : Base(ham) {}

        protected:
            virtual void singles(conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void singles(const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void doubles(conn::FrmOnv &conn, const std::function<void()> &fn);

            virtual void foreach(const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn);

        public:
            void foreach(const FrmOnv &mbf, const fn_c_t<FrmOnv> &body_fn) override;

            void foreach(const FrmBosOnv &mbf, const fn_c_t<FrmBosOnv> &body_fn) override;

            void foreach(const BosOnv &mbf, const fn_c_t<BosOnv> &body_fn) override {}

        };

        /**
         * make use of spin and point-group symmetry in the loop over connections
         */
        struct SpinSym : Fermion {
            using Fermion::foreach;

            explicit SpinSym(const Hamiltonian &ham) : Fermion(ham) {}

        protected:
            void singles(conn::FrmOnv &conn, const std::function<void()> &fn) override;

            void doubles(conn::FrmOnv &conn, const std::function<void()> &fn) override;
        };

        struct Hubbard1D : Fermion {
            using Fermion::foreach;
            const bool m_pbc;

            explicit Hubbard1D(const Hamiltonian &ham);

        protected:
            void singles(const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) override;

        public:

            void foreach(const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) override;
        };
    }

    static std::unique_ptr<Base> make(const Hamiltonian& ham) {
        return std::unique_ptr<Base>(new frm::Fermion(ham));
    }
}

#endif //M7_FOREACHCONNECTION_H
