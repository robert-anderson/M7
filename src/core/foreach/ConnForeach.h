//
// Created by anderson on 17/02/2022.
//

#ifndef M7_CONNFOREACH_H
#define M7_CONNFOREACH_H

#include <utility>
#include <src/core/caches/CachedOrbs.h>

#include "src/core/basis/Suites.h"
#include "src/core/connection/Connections.h"
#include "ForeachVirtual.h"

#if 0
/**
 * Foreach-type iterators for inter-MBF connections of a certain excitation signature
 */
namespace conn_foreach {

    using namespace foreach_virtual::ctnd;

    struct ConnForeach {
        const BasisDims m_bd;
        CachedOrbs m_work_orbs;

        explicit ConnForeach(BasisDims bd) : m_bd(bd){}

        ConnForeach(const ConnForeach &other) : ConnForeach(other.m_bd){}
        virtual ~ConnForeach(){}

        virtual void throwing_loop() = 0;

        /**
         * iterates over all values and iiters
         */
        void loop(){
            try {throwing_loop();}
            catch (const ExitLoop&){}
        }

        /**
         * @return
         *  current value of the iteration counter
         */
        virtual size_t iiter() const = 0;

        /**
         * @return
         *  number of calls to body if no early termination
         */
        virtual size_t niter() const = 0;
    };


    template<size_t mbf_ind>
    class Base : public ConnForeach {
    protected:
        typedef conn::mbf_t<mbf_ind> conn_t;
        conn_t m_conn_internal;
        conn_t *m_conn;

        typedef std::function<void(const conn_t&, size_t)> body_fn_t;
        body_fn_t m_body_fn;
    public:
        Base(BasisDims bd, body_fn_t body_fn = {}, conn_t *conn = nullptr) :
                ConnForeach(bd), m_conn_internal(bd), m_conn(conn ? conn : &m_conn_internal),
                m_body_fn(std::move(body_fn)) {}

        Base(const Base &other, conn_t *conn = nullptr) : Base(other.m_bd, other.m_body_fn, conn) {}

        virtual void body() {
            if (m_body_fn) m_body_fn(*m_conn, iiter());
        }

        const conn_t& conn() const {
            return *m_conn;
        }
    };


    namespace frm {
        class Base : public conn_foreach::Base<defs::Frm> {
            const field::FrmOnv& m_src;
        public:
            Base(const field::FrmOnv& src, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                conn_foreach::Base<defs::Frm>({src.nsite(), 0ul}, std::move(body_fn), conn),
                m_src(src){}

            void frm_throwing_loop(const field::FrmOnv &src) override {
                m_work_orbs.clear();
            }
            void bos_throwing_loop(const field::BosOnv &src) override {}
            void frmbos_throwing_loop(const field::FrmBosOnv &src) override {
                frm_throwing_loop(src.m_frm);
            }
        };

        template<size_t nexcit>
        class General : public Base {
            const defs::inds& m_occ;
            const defs::inds& m_vac;

            struct Foreach : Ordered<nexcit> {
                General& m_context;
                Foreach(General& context):
                    Ordered<nexcit>((context.m_work_orbs.nbit), m_context(context){}

                void body() override {
                    //m_context.m_conn->set()
                }
            };
            Foreach m_foreach;

        public:
            General(size_t nsite, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                Base(nsite, std::move(body_fn), conn),
                m_occ(m_work_orbs.occ(mbf)),
                m_foreach(*this, 2*nsite){}

            void frm_throwing_loop(const field::FrmOnv &src) override {
                Base::frm_throwing_loop(src);
                m_foreach.loop();
            }

            size_t iiter() const override {
                return 0;
            }

            size_t niter() const override {
                return 0;
            }
        };
    }

    namespace bos {
        class Base : public conn_foreach::Base<defs::Bos> {
        public:
            Base(size_t nmode, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                conn_foreach::Base<defs::Bos>({0ul, nmode}, std::move(body_fn), conn){}

            void frm_throwing_loop(const field::FrmOnv &src) override {}
            void frmbos_throwing_loop(const field::FrmBosOnv &src) override {
                bos_throwing_loop(src.m_bos);
            }
        };
    }

    namespace frm_bos {

        class Base : public conn_foreach::Base<defs::FrmBos> {
        public:
            Base(BasisDims bd, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                conn_foreach::Base<defs::FrmBos>(bd, std::move(body_fn), conn){}

            void frm_throwing_loop(const field::FrmOnv &src) override {}
            void bos_throwing_loop(const field::BosOnv &src) override {}
        };
    }


#if 0
    namespace frm {}
    namespace bos {}
    namespace frm_bos {}


}
class ConnForeach {
#endif

};


#endif //M7_CONNFOREACH_H
#endif //M7_CONNFOREACH_H
