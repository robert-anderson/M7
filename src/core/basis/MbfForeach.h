//
// Created by anderson on 2/9/22.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H


#include <src/core/field/Fields.h>
#include <src/core/util/ForeachVirtual.h>

namespace mbf_foreach {
    using namespace foreach_virtual::rtnd;

    struct FrmOnv {
        field::FrmOnv &m_mbf;
        FrmOnv(field::FrmOnv &mbf) : m_mbf(mbf) {}
        virtual void body() = 0;
        virtual void loop() = 0;
    };

//    struct FrmOnvSpinConserve : FrmOnv {
//        /**
//         * iterator over a single spin channel
//         */
//        struct Foreach : Ordered<> {
//            Foreach(size_t nsite, size_t nelec): Ordered<>(nsite, )
//        };
//    };

}

#if 0
namespace mbf_foreach {


    struct FrmOnv {
        field::FrmOnv &m_frmonv;
        FrmOnv(field::FrmOnv &mbf) : m_frmonv(mbf) {}
        virtual void mbf_body(const field::FrmOnv& mbf) = 0;
    };
    struct FrmBosOnv {
        field::FrmBosOnv &m_frmbosonv;
        FrmBosOnv(field::FrmBosOnv &mbf) : m_frmbosonv(mbf) {}
        virtual void mbf_body(const field::FrmBosOnv& mbf) = 0;
    };

    struct BosOnv {
        field::BosOnv &m_bosonv;
        BosOnv(field::BosOnv &mbf) : m_bosonv(mbf) {}
        virtual void mbf_body(const field::BosOnv& mbf) = 0;
    };

    struct FrmOnvGeneral : FrmOnv, Ordered<true, true> {

        FrmOnvGeneral(field::FrmOnv &mbf, size_t nelec): FrmOnv(mbf), Ordered<true, true>(mbf.m_nspinorb, nelec){}

        void body() override {
            m_frmonv = inds();
            mbf_body(m_frmonv);
        }
    };

    struct BosOnvGeneral : BosOnv, Unrestricted {

        BosOnvGeneral(field::BosOnv &mbf, size_t nboson_max):
            BosOnv(mbf), Unrestricted(mbf.m_nmode, nboson_max+1){}

        void body() override {
            m_bosonv = inds();
            mbf_body(m_bosonv);
        }
    };


    struct FrmBosOnvGeneral : FrmBosOnv {
        field::FrmBosOnv& m_mbf;
        const size_t m_nelec;
        const size_t m_nboson_max;

        struct Inner : BosOnvGeneral {
            FrmBosOnvGeneral& m_context;
            Inner(FrmBosOnvGeneral& context):
                BosOnvGeneral(context.m_mbf.m_bos, context.m_nboson_max), m_context(context){}

            void mbf_body(const field::BosOnv &mbf) override {
                m_context.mbf_body(m_context.m_mbf);
            }
        };

        struct Outer : FrmOnvGeneral {
            FrmBosOnvGeneral& m_context;
            Outer(FrmBosOnvGeneral& context):
                FrmOnvGeneral(context.m_mbf.m_frm, context.m_nelec), m_context(context){}

            void mbf_body(const field::FrmOnv &mbf) override {
                Inner inner(m_context);
                inner.loop();
            }
        };

        FrmBosOnvGeneral(field::FrmBosOnv& mbf, size_t nelec, size_t nboson_max):
            FrmBosOnv(mbf), m_mbf(mbf), m_nelec(nelec), m_nboson_max(nboson_max){}

//        void mbf_body(const field::FrmBosOnv &mbf) override {
//            std::cout << m_frmonv << std::endl;
//        }

        void loop() {
            Outer(*this).loop();
        }
    };


    template<typename frm_t, typename bos_t>
    struct FrmBosProduct : FrmBosOnv {
        field::FrmBosOnv &m_mbf;
        const size_t m_nelec;
        const size_t m_nboson_max;

        struct Inner : BosOnvGeneral {
            FrmBosOnvGeneral &m_context;

            Inner(FrmBosOnvGeneral &context) :
                    BosOnvGeneral(context.m_mbf.m_bos, context.m_nboson_max), m_context(context) {}

            void mbf_body(const field::BosOnv &mbf) override {
                m_context.mbf_body(m_context.m_mbf);
            }
        };

        struct Outer : FrmOnvGeneral {
            FrmBosOnvGeneral &m_context;

            Outer(FrmBosOnvGeneral &context) :
                    FrmOnvGeneral(context.m_mbf.m_frm, context.m_nelec), m_context(context) {}

            void mbf_body(const field::FrmOnv &mbf) override {
                Inner inner(m_context);
                inner.loop();
            }
        };

        FrmBosProduct(field::FrmBosOnv &mbf, size_t nelec, size_t nboson_max) :

        FrmBosOnv (mbf), m_mbf(mbf), m_nelec(nelec), m_nboson_max(nboson_max) {}

        void loop() {
            Outer(*this).loop();
        }
    };
}


#endif //M7_MBFFOREACH_H
#endif //M7_MBFFOREACH_H
