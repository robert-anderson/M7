//
// Created by anderson on 2/9/22.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H


#include <src/core/field/Fields.h>
#include <src/core/util/ForeachVirtual.h>
#include <src/core/table/BufferedFields.h>

namespace mbf_foreach {
    using namespace foreach_virtual::rtnd;

    struct FrmOnv {
        buffered::FrmOnv m_mbf;
        FrmOnv(size_t nsite) : m_mbf(nsite) {}
        virtual void body() = 0;
        virtual void loop() = 0;
        virtual size_t iiter() = 0;
    };

    struct FrmOnvGeneral : FrmOnv {
        struct Foreach : Ordered<> {
            FrmOnvGeneral& m_context;
            Foreach(FrmOnvGeneral& context, size_t nelec):
                    Ordered<>(context.m_mbf.m_nspinorb, nelec), m_context(context){}

            void body() override {
                field::NdBitset<size_t, 2>& mbf = m_context.m_mbf;
                mbf.zero();
                mbf = value();
                m_context.body();
            }
        };
        Foreach m_foreach;
        FrmOnvGeneral(size_t nsite, size_t nelec): FrmOnv(nsite), m_foreach(*this, nelec){}

        void loop() override {
            m_foreach.loop();
        }

        size_t iiter() override {
            return m_foreach.iiter();
        }
    };

    struct FrmOnvSpins : FrmOnv {
        struct Foreach : Ordered<> {
            FrmOnvSpins& m_context;
            Foreach(FrmOnvSpins& context, int ms2):
                    Ordered<>(context.m_mbf.m_nsite, (context.m_mbf.m_nsite+ms2)/2), m_context(context){}

            void body() override {
                m_context.m_mbf.set_spins(value());
                m_context.body();
            }
        };
        Foreach m_foreach;
        FrmOnvSpins(size_t nsite, int ms2): FrmOnv(nsite), m_foreach(*this, ms2){}

        void loop() override {
            m_foreach.loop();
        }

        size_t iiter() override {
            return m_foreach.iiter();
        }
    };

    struct FrmOnvSzConserve : FrmOnv {
        /**
         * iterator over a single spin channel
         */
        struct Foreach : Ordered<> {
            FrmOnvSzConserve& m_context;
            Foreach(FrmOnvSzConserve& context, size_t nelec):
                Ordered<>(context.m_mbf.nsite(), nelec), m_context(context){}
        };
        /**
         * iterator over inner spin channel
         */
        struct BetaForeach : Foreach {
            BetaForeach(FrmOnvSzConserve& context, size_t nelec): Foreach(context, nelec){}
            void body() override {
                m_context.m_mbf.put_spin_channel(1, false);
                m_context.m_mbf.set(m_context.m_mbf.m_nsite, value());
                m_context.body();
            }
        };
        /**
         * iterator over outer spin channel
         */
        struct AlphaForeach : Foreach {
            AlphaForeach(FrmOnvSzConserve& context, size_t nelec): Foreach(context, nelec){}
            void body() override {
                m_context.m_mbf.put_spin_channel(0, false);
                m_context.m_mbf.set(0, value());
                m_context.m_beta_foreach.loop();
            }
        };

        AlphaForeach m_alpha_foreach;
        BetaForeach m_beta_foreach;
        FrmOnvSzConserve(size_t nsite, size_t nelec): FrmOnv(nsite),
                                                      m_alpha_foreach(*this, nelec), m_beta_foreach(*this, nelec){}

        void loop() override {
            m_alpha_foreach.loop();
        }

        size_t iiter() override {
            return m_alpha_foreach.iiter()*m_beta_foreach.m_niter + m_beta_foreach.iiter();
        }
    };

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
