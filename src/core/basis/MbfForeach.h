//
// Created by anderson on 2/9/22.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H


#include <src/core/field/Fields.h>
#include <src/core/util/ForeachVirtual.h>
#include <src/core/table/BufferedFields.h>

namespace mbf_foreach {
    /**
     * all MBFs have runtime-determined dimensions
     */
    using namespace foreach_virtual::rtnd;

    struct Base {
        const BasisDims m_bd;

        Base(BasisDims bd) : m_bd(bd) {}

        /**
         * function to inject into structured loop, referencing value and iiter methods just as in the primitive
         * foreach iterators
         */
        virtual void body() = 0;

        /**
         * iterates over all values and iiters
         */
        virtual void loop() = 0;

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

    namespace frm {

        class Base : public mbf_foreach::Base {
            buffered::FrmOnv m_mbf_internal;
        protected:
            field::FrmOnv *m_mbf;
        public:
            Base(size_t nsite, field::FrmOnv *mbf = nullptr) :
                    mbf_foreach::Base({nsite, 0ul}),
                    m_mbf_internal(nsite), m_mbf(mbf ? mbf : &m_mbf_internal) {}

            Base(const Base &other, field::FrmOnv *mbf = nullptr) : Base(other.m_mbf->m_nsite, mbf) {}

            const field::FrmOnv &value() const {
                return *m_mbf;
            }
        };

        class General : public Base {
            struct Foreach : Ordered<> {
                General &m_context;

                Foreach(General &context, size_t nelec) :
                        Ordered<>(context.m_mbf->m_nspinorb, nelec), m_context(context) {}

                void body() override {
                    m_context.m_mbf->zero();
                    *m_context.m_mbf = value();
                    m_context.body();
                }
            };

            Foreach m_foreach;

        public:
            General(size_t nsite, size_t nelec, field::FrmOnv *mbf = nullptr) :
                    Base(nsite, mbf), m_foreach(*this, nelec) {}

            General(const General &other, field::FrmOnv *mbf = nullptr) :
                    General(other.m_bd.m_nsite, other.m_foreach.m_nind, mbf) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() const override {
                return m_foreach.iiter();
            }

            size_t niter() const override {
                return m_foreach.m_niter;
            }
        };

        class Spins : public Base {
            struct Foreach : Ordered<> {
                Spins &m_context;

                Foreach(Spins &context, int ms2) :
                        Ordered<>(context.m_mbf->m_nsite, (context.m_mbf->m_nsite + ms2) / 2), m_context(context) {}

                void body() override {
                    m_context.m_mbf->set_spins(value());
                    m_context.body();
                }

            };

            Foreach m_foreach;

            int ms2() const {
                return 2 * m_foreach.m_nind - m_mbf->m_nsite;
            }

        public:
            Spins(size_t nsite, int ms2, field::FrmOnv *mbf = nullptr) : Base(nsite, mbf), m_foreach(*this, ms2) {}

            Spins(const Spins &other, field::FrmOnv *mbf = nullptr) : Spins(other.m_bd.m_nsite, other.ms2(), mbf) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() const override {
                return m_foreach.iiter();
            }

            size_t niter() const override {
                return m_foreach.m_niter;
            }
        };

        using namespace ci_utils;

        class Ms2Conserve : public Base {
            /**
             * iterator over a single spin channel
             */
            struct Foreach : Ordered<> {
                Ms2Conserve &m_context;

                Foreach(Ms2Conserve &context, size_t nelec) :
                        Ordered<>(context.m_mbf->m_nsite, nelec), m_context(context) {}
            };

            /**
             * iterator over inner spin channel
             */
            struct BetaForeach : Foreach {
                BetaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body() override {
                    m_context.m_mbf->put_spin_channel(1, false);
                    m_context.m_mbf->set(m_context.m_mbf->m_nsite, value());
                    m_context.body();
                }
            };

            /**
             * iterator over outer spin channel
             */
            struct AlphaForeach : Foreach {
                AlphaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body() override {
                    m_context.m_mbf->put_spin_channel(0, false);
                    m_context.m_mbf->set(0, value());
                    m_context.m_beta_foreach.loop();
                }
            };

            AlphaForeach m_alpha_foreach;
            BetaForeach m_beta_foreach;


            size_t nelec() const {
                return m_alpha_foreach.m_nind + m_beta_foreach.m_nind;
            }

            int ms2() const {
                return m_alpha_foreach.m_nind - m_beta_foreach.m_nind;
            }

        public:
            Ms2Conserve(size_t nsite, size_t nelec, int ms2, field::FrmOnv *mbf = nullptr) :
                    Base(nsite, mbf),
                    m_alpha_foreach(*this, nalpha(nelec, ms2)),
                    m_beta_foreach(*this, nbeta(nelec, ms2)) {}

            Ms2Conserve(const Ms2Conserve &other, field::FrmOnv *mbf = nullptr) :
                    Ms2Conserve(other.m_bd.m_nsite, other.nelec(), other.ms2(), mbf) {}

            void loop() override {
                m_alpha_foreach.loop();
            }

            size_t iiter() const override {
                return m_alpha_foreach.iiter() * m_beta_foreach.m_niter + m_beta_foreach.iiter();
            }

            size_t niter() const override {
                return m_alpha_foreach.m_niter * m_beta_foreach.m_niter;
            }
        };
    }

    namespace bos {
        class Base : public mbf_foreach::Base {
            buffered::BosOnv m_mbf_internal;
        protected:
            field::BosOnv *m_mbf;
        public:
            Base(size_t nmode, field::BosOnv *mbf = nullptr) :
                    mbf_foreach::Base({0, nmode}),
                    m_mbf_internal(nmode), m_mbf(mbf ? mbf : &m_mbf_internal) {}

            Base(const Base &other, field::BosOnv *mbf = nullptr) : Base(other.m_bd.m_nmode, mbf) {}

            const field::BosOnv &value() const {
                return *m_mbf;
            }
        };

        class General : public Base {
            struct Foreach : Unrestricted {
                General &m_context;

                Foreach(General &context, size_t nboson_max) :
                        Unrestricted(context.m_mbf->m_nmode, nboson_max + 1), m_context(context) {}

                void body() override {
                    *m_context.m_mbf = value();
                    m_context.body();
                }

                size_t nboson_max() const {
                    if (m_shape.empty()) return ~0ul;
                    return m_shape[0] - 1;
                }
            };

            Foreach m_foreach;
        public:
            General(size_t nmode, size_t nboson_max, field::BosOnv *mbf = nullptr) :
                    Base(nmode, mbf), m_foreach(*this, nboson_max) {}

            General(const General &other, field::BosOnv *mbf = nullptr) :
                    General(other.m_bd.m_nmode, other.m_foreach.nboson_max(), mbf) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() const override {
                return m_foreach.iiter();
            }

            size_t niter() const override {
                return m_foreach.m_niter;
            }
        };
    }

    namespace frm_bos {

        class Base : public mbf_foreach::Base {
            buffered::FrmBosOnv m_mbf_internal;
        protected:
            field::FrmBosOnv *m_mbf;
        public:
            Base(BasisDims bd, field::FrmBosOnv *mbf = nullptr) :
                    mbf_foreach::Base(bd), m_mbf_internal(bd), m_mbf(mbf ? mbf : &m_mbf_internal) {}

            Base(const Base &other, field::FrmBosOnv *mbf = nullptr) : Base(other.m_bd, mbf) {}

            const field::FrmBosOnv &value() const {
                return *m_mbf;
            }
        };

        template<typename frm_foreach_t, typename bos_foreach_t>
        class Product : public Base {
            static_assert(std::is_base_of<mbf_foreach::frm::Base, frm_foreach_t>::value,
                          "template arg must be derived from the base type of fermion foreach iterators");
            static_assert(std::is_base_of<mbf_foreach::bos::Base, bos_foreach_t>::value,
                          "template arg must be derived from the base type of boson foreach iterators");

            struct BosForeach : public bos_foreach_t {
                Product &m_context;

                BosForeach(Product &context, const bos_foreach_t &other) :
                        bos_foreach_t(other, &context.m_mbf->m_bos), m_context(context) {
                    REQUIRE_EQ(bos::Base::m_mbf, &m_context.m_mbf->m_bos, "wrong pointer");
                }

                void body() override {
                    bos_foreach_t::body();
                    m_context.body();
                }
            };

            struct FrmForeach : public frm_foreach_t {
                Product &m_context;

                FrmForeach(Product &context, const frm_foreach_t &other) :
                        frm_foreach_t(other, &context.m_mbf->m_frm), m_context(context) {
                    REQUIRE_EQ(frm::Base::m_mbf, &m_context.m_mbf->m_frm, "wrong pointer");
                }

                void body() override {
                    frm_foreach_t::body();
                    m_context.m_bos_foreach.loop();
                }
            };

        public:
            FrmForeach m_frm_foreach;
            BosForeach m_bos_foreach;

            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach, field::FrmBosOnv *mbf=nullptr):
                Base({static_cast<const mbf_foreach::frm::Base &>(frm_foreach).m_bd.m_nsite,
                        static_cast<const mbf_foreach::bos::Base &>(bos_foreach).m_bd.m_nmode}, mbf),
                    m_frm_foreach(*this, frm_foreach), m_bos_foreach(*this, bos_foreach) {}

            void loop() override {
                m_frm_foreach.loop();
            }

            size_t iiter() const override {
                auto iiter_frm = static_cast<const mbf_foreach::Base &>(m_frm_foreach).iiter();
                auto iiter_bos = static_cast<const mbf_foreach::Base &>(m_bos_foreach).iiter();
                auto niter_bos = static_cast<const mbf_foreach::Base &>(m_bos_foreach).niter();
                return iiter_frm * niter_bos + iiter_bos;
            }

            size_t niter() const override {
                auto niter_frm = static_cast<const mbf_foreach::Base &>(m_frm_foreach).niter();
                auto niter_bos = static_cast<const mbf_foreach::Base &>(m_bos_foreach).niter();
                return niter_frm * niter_bos;
            }
        };
    }
}

#endif //M7_MBFFOREACH_H
