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

    template<typename mbf_t>
    using fn_cd_t = std::function<void(const conn::from_field_t<mbf_t>&, const mbf_t&)>;
    template<typename mbf_t>
    using fn_ch_t = std::function<void(const conn::from_field_t<mbf_t>&, defs::ham_t)>;
    template<typename mbf_t>
    using fn_cdh_t = std::function<void(const conn::from_field_t<mbf_t>&, const mbf_t&, defs::ham_t)>;
    template<typename mbf_t>
    using fn_d_t = std::function<void(const mbf_t&)>;
    using fn_h_t = std::function<void(defs::ham_t)>;
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

        /*
         * these virtual methods define what to do with each connection of the given many-body basis function
         */
        virtual void foreach(const FrmOnv& mbf, const fn_c_t<FrmOnv>& body_fn) = 0;
        virtual void foreach(const FrmBosOnv& mbf, const fn_c_t<FrmBosOnv>& body_fn) = 0;
        virtual void foreach(const BosOnv& mbf, const fn_c_t<BosOnv>& body_fn) = 0;

        /*
         * The following definitions adapt the body functions of various types for use with the above virtual methods.
         *
         * For now, these are kept private since the user can just call one of the foreach overloads and the adapt
         * method will be automatically called before delegating to the above fn_c_t foreach definition. There is a
         * case for encouraging the use of the adapt method in the calling scope of foreach, since the same wrapper
         * std::function object is otherwise being constructed on every call to the foreach overload. However, this
         * could create more object lifetime problems if used incorrectly, and in any case, the creation of said
         * wrappers should be a small overhead compared to the time taken to execute the foreach method itself.
         */
    private:
        /**
         * adapts "c" closure to "c" closure which optionally checks for nonzero H matrix element before calling body_fn
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connection and connected MBF args
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_c_t<mbf_t>& body_fn, bool nonzero_h_only){
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                if (nonzero_h_only && consts::float_is_zero(get_element(mbf, conn))) return;
                body_fn(conn);
            };
        }
        /**
         * adapts "cd" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connection and connected MBF args
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_cd_t<mbf_t>& body_fn, bool nonzero_h_only){
            auto& dst_mbf = m_mbfs[mbf];
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                if (nonzero_h_only && consts::float_is_zero(get_element(mbf, conn))) return;
                conn.apply(mbf, dst_mbf);
                body_fn(conn, dst_mbf);
            };
        }

        /**
         * adapts "ch" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connection and H matrix element args
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_ch_t<mbf_t>& body_fn, bool nonzero_h_only){
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                auto helem = get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                body_fn(conn, helem);
            };
        }

        /**
         * adapts "cdh" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connection, connected MBF args, and H matrix element as args
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_cdh_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto &dst_mbf = m_mbfs[mbf];
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                auto helem = get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                conn.apply(mbf, dst_mbf);
                body_fn(conn, dst_mbf, helem);
            };
        }

        /**
         * adapts "d" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connected MBF as its only argument
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_d_t<mbf_t>& body_fn, bool nonzero_h_only){
            auto& dst_mbf = m_mbfs[mbf];
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                if (nonzero_h_only && consts::float_is_zero(get_element(mbf, conn))) return;
                conn.apply(mbf, dst_mbf);
                body_fn(dst_mbf);
            };
        }

        /**
         * adapts "h" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts H matrix element as its only argument
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_h_t& body_fn){
            return [&](const conn::from_field_t<mbf_t> &conn) {
                auto helem = get_element(mbf, conn);
                if (consts::float_is_zero(helem)) return;
                body_fn(helem);
            };
        }

        /**
         * adapts "dh" closure to "c" closure
         * @tparam mbf_t
         *  many-body basis function type
         * @param mbf
         *  source many-body basis function from which to enumerate all connections and call the given function
         * @param body_fn
         *  void function which accepts connection, connected MBF args, and H matrix element as args
         * @param nonzero_h_only
         *  if true, body_fn is only called when the loop generates a connection with non-zero H matrix element
         */
        template<typename mbf_t>
        fn_c_t<mbf_t> adapt(const mbf_t& mbf, const fn_dh_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto &dst_mbf = m_mbfs[mbf];
            return [&, nonzero_h_only](const conn::from_field_t<mbf_t> &conn) {
                auto helem = get_element(mbf, conn);
                if (nonzero_h_only && consts::float_is_zero(helem)) return;
                conn.apply(mbf, dst_mbf);
                body_fn(dst_mbf, helem);
            };
        }

    public:
        /*
         * The following definitions adapt the body functions of various types for use with the above "fn_c_t" foreach
         * virtual methods
         */

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_c_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_cd_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_ch_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_cdh_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_d_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_h_t& body_fn) {
            auto fn = adapt(mbf, body_fn);
            this->foreach(mbf, fn);
        }

        template<typename mbf_t>
        void foreach(const mbf_t& mbf, const fn_dh_t<mbf_t>& body_fn, bool nonzero_h_only) {
            auto fn = adapt(mbf, body_fn, nonzero_h_only);
            this->foreach(mbf, fn);
        }


        defs::ham_t get_element(const FrmOnv& mbf, const conn::FrmOnv& conn);
        defs::ham_t get_element(const FrmBosOnv& mbf, const conn::FrmBosOnv& conn);
        defs::ham_t get_element(const BosOnv& mbf, const conn::BosOnv& conn);
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
