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

    struct Base {
        virtual void body() = 0;

        virtual void loop() = 0;

        virtual size_t iiter() = 0;
    };

    namespace frm {

        struct Base : mbf_foreach::Base {
            buffered::FrmOnv m_mbf;
            Base(size_t nsite) : m_mbf(nsite) {}
        };

        struct General : Base {
            struct Foreach : Ordered<> {
                General &m_context;

                Foreach(General &context, size_t nelec) :
                        Ordered<>(context.m_mbf.m_nspinorb, nelec), m_context(context) {}

                void body() override {
                    field::NdBitset<size_t, 2> &mbf = m_context.m_mbf;
                    mbf.zero();
                    mbf = value();
                    m_context.body();
                }
            };

            Foreach m_foreach;

            General(size_t nsite, size_t nelec) : Base(nsite), m_foreach(*this, nelec) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
            }
        };

        struct Spins : Base {
            struct Foreach : Ordered<> {
                Spins &m_context;

                Foreach(Spins &context, int ms2) :
                        Ordered<>(context.m_mbf.m_nsite, (context.m_mbf.m_nsite + ms2) / 2), m_context(context) {}

                void body() override {
                    m_context.m_mbf.set_spins(value());
                    m_context.body();
                }
            };

            Foreach m_foreach;

            Spins(size_t nsite, int ms2) : Base(nsite), m_foreach(*this, ms2) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
            }
        };

        struct Ms2Conserve : Base {
            /**
             * iterator over a single spin channel
             */
            struct Foreach : Ordered<> {
                Ms2Conserve &m_context;

                Foreach(Ms2Conserve &context, size_t nelec) :
                        Ordered<>(context.m_mbf.nsite(), nelec), m_context(context) {}
            };

            /**
             * iterator over inner spin channel
             */
            struct BetaForeach : Foreach {
                BetaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

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
                AlphaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body() override {
                    m_context.m_mbf.put_spin_channel(0, false);
                    m_context.m_mbf.set(0, value());
                    m_context.m_beta_foreach.loop();
                }
            };

            AlphaForeach m_alpha_foreach;
            BetaForeach m_beta_foreach;

            Ms2Conserve(size_t nsite, size_t nelec, int ms2) :
                    Base(nsite),
                    m_alpha_foreach(*this, ci_utils::nalpha(nelec, ms2)),
                    m_beta_foreach(*this, ci_utils::nbeta(nelec, ms2)) {}

            void loop() override {
                m_alpha_foreach.loop();
            }

            size_t iiter() override {
                return m_alpha_foreach.iiter() * m_beta_foreach.m_niter + m_beta_foreach.iiter();
            }
        };
    }

    namespace bos {
        struct Base : mbf_foreach::Base {
            buffered::BosOnv m_mbf;
            Base(size_t nmode) : m_mbf(nmode) {}
        };

        struct General : Base {
            struct Foreach : Unrestricted {
                General &m_context;

                Foreach(General &context, size_t nboson_max) :
                        Unrestricted(context.m_mbf.m_nmode, nboson_max), m_context(context) {}

                void body() override {
                    m_context.m_mbf = value();
                }
            };

            Foreach m_foreach;
            General(size_t nmode, size_t nboson_max): Base(nmode), m_foreach(*this, nboson_max){}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
            }
        };
    }

    namespace frmbos {

        struct Base : mbf_foreach::Base {
            buffered::FrmBosOnv m_mbf;
            Base(BasisDims bd) : m_mbf(bd) {}
        };

    }

}

#endif //M7_MBFFOREACH_H
