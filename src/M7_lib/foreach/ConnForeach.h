//
// Created by anderson on 17/02/2022.
//

#ifndef M7_CONNFOREACH_H
#define M7_CONNFOREACH_H

#include <utility>
#include <array>

#include <M7_lib/basis/Suites.h>
#include <M7_lib/connection/Connections.h>

#include "ForeachVirtual.h"


namespace conn_foreach {
    using namespace foreach_virtual::ctnd;

    struct ConnForeach {
        const size_t m_exsig;
        ConnForeach(size_t exsig): m_exsig(exsig){}
        virtual void frm_throwing_loop(const FrmOnvField &src) = 0;
        virtual void bos_throwing_loop(const BosOnvField &src) = 0;
        virtual void frmbos_throwing_loop(const FrmBosOnvField &src) = 0;
        void throwing_loop(const FrmOnvField& src){
            frm_throwing_loop(src);
        }
        void throwing_loop(const BosOnvField& src){
            bos_throwing_loop(src);
        }
        void throwing_loop(const FrmBosOnvField& src){
            frmbos_throwing_loop(src);
        }

        virtual size_t frm_niter(const FrmOnvField &src) = 0;
        virtual size_t bos_niter(const BosOnvField &src) = 0;
        virtual size_t frmbos_niter(const FrmBosOnvField &src) = 0;

        size_t niter(const FrmOnvField &src) {
            return frm_niter(src);
        }
        size_t niter(const BosOnvField &src) {
            return bos_niter(src);
        }
        size_t niter(const FrmBosOnvField &src) {
            return frmbos_niter(src);
        }

        template<class mbf_t>
        void loop(const mbf_t& src){
            try {throwing_loop(src);}
            catch (const ExitLoop&){}
        }
    };


    template<size_t mbf_ind>
    class Base : public ConnForeach {
    protected:
        typedef conn::mbf_t<mbf_ind> conn_t;
        conn_t m_conn_internal;
        conn_t *m_conn;

        typedef std::function<void(const conn_t&)> body_fn_t;
        body_fn_t m_body_fn;
    public:
        Base(size_t exsig, BasisData bd, body_fn_t body_fn = {}, conn_t *conn = nullptr) :
                ConnForeach(exsig), m_conn_internal(bd),
                m_conn(conn ? conn : &m_conn_internal), m_body_fn(std::move(body_fn)) {}

        Base(const Base &other, conn_t *conn = nullptr) : Base(other.m_bd, other.m_body_fn, conn) {}

        virtual void body() {
            if (m_body_fn) m_body_fn(*m_conn);
        }
    };

    namespace frm {
        class Base : public conn_foreach::Base<defs::Frm> {
        public:
            Base(size_t exsig, size_t nsite, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                    conn_foreach::Base<defs::Frm>(exsig, {nsite, 0ul}, std::move(body_fn), conn){}

            void bos_throwing_loop(const field::BosOnv &src) override {}
            void frmbos_throwing_loop(const field::FrmBosOnv &src) override {
                frm_throwing_loop(src.m_frm);
            }

            size_t bos_niter(const BosOnvField &src) override {
                return 0;
            }

            size_t frmbos_niter(const FrmBosOnvField &src) override {
                return frm_niter(src.m_frm);
            }
        };

        template<size_t nexcit>
        class General : public Base {

            struct VacForeach : Ordered<nexcit> {
                General& m_context;
                const defs::inds& m_vacs;
                VacForeach(General& context, const defs::inds& vacs):
                        Ordered<nexcit>(vacs.size()), m_context(context), m_vacs(vacs){}

                void body(const inds_t<nexcit> &value, size_t iiter) override {
                    m_context.m_conn->m_cre.clear();
                    for (size_t i=0ul; i<nexcit; ++i) m_context.m_conn->m_cre.add(m_vacs[value[i]]);
                    m_context.body();
                }
            };

            struct OccForeach : Ordered<nexcit> {
                General& m_context;
                const defs::inds& m_occs;
                const defs::inds& m_vacs;
                OccForeach(General& context, const defs::inds& occs, const defs::inds& vacs):
                        Ordered<nexcit>(occs.size()), m_context(context), m_occs(occs), m_vacs(vacs){}

                void body(const inds_t<nexcit> &value, size_t iiter) override {
                    m_context.m_conn->clear();
                    for (size_t i=0ul; i<nexcit; ++i) m_context.m_conn->m_ann.add(m_occs[value[i]]);
                    VacForeach(m_context, m_vacs).loop();
                }
            };

        public:
            General(size_t nsite, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                    Base(exsig_utils::encode(nexcit, nexcit, 0, 0), nsite, std::move(body_fn), conn){}

            void frm_throwing_loop(const field::FrmOnv &src) override {
                src.m_decoded.clear(); // reset all caches
                auto &occs = src.m_decoded.m_simple_occs.get();
                auto &vacs = src.m_decoded.m_simple_vacs.get();
                if (occs.empty() || vacs.empty()) return;
                OccForeach(*this, occs, vacs).loop();
            }

            size_t frm_niter(const FrmOnvField &src) override {
                auto &occs = src.m_decoded.m_simple_occs.get();
                auto &vacs = src.m_decoded.m_simple_vacs.get();
                auto nocc_comb = integer_utils::combinatorial(occs.size(), nexcit);
                auto nvac_comb = integer_utils::combinatorial(vacs.size(), nexcit);
                return nocc_comb*nvac_comb;
            }
        };
    }

    namespace bos {

    }
}

#if 0
/**
 * Foreach-type iterators for inter-MBF connections of a certain excitation signature
 */
namespace conn_foreach {

    using namespace foreach_virtual::ctnd;

    struct ConnForeach {
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
        Base(BasisData bd, body_fn_t body_fn = {}, conn_t *conn = nullptr) :
                m_conn_internal(bd), m_conn(conn ? conn : &m_conn_internal),
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

#if 0
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
            Base(BasisData bd, body_fn_t body_fn = {}, conn_t *conn = nullptr):
                conn_foreach::Base<defs::FrmBos>(bd, std::move(body_fn), conn){}

            void frm_throwing_loop(const field::FrmOnv &src) override {}
            void bos_throwing_loop(const field::BosOnv &src) override {}
        };
    }


    namespace frm {}
    namespace bos {}
    namespace frm_bos {}


}
class ConnForeach {
#endif

}


#endif //M7_CONNFOREACH_H
#endif //M7_CONNFOREACH_H
