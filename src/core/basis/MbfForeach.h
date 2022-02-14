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
        virtual size_t iiter() = 0;
    };

    namespace frm {

        class Base : public mbf_foreach::Base {
            buffered::FrmOnv m_mbf_internal;
        protected:
            field::FrmOnv &m_mbf;
        public:
            Base(size_t nsite) : mbf_foreach::Base({nsite, 0ul}), m_mbf_internal(nsite), m_mbf(m_mbf_internal) {}

            Base(field::FrmOnv &mbf) : mbf_foreach::Base({mbf.m_nsite, 0ul}), m_mbf_internal(mbf), m_mbf(mbf) {}

            /**
             * the m_mbf reference-switching copy ctor is trivial here, but its use becomes clear in the FrmBos products
             * @param mbf
             *  new reference to be used
             * @param other
             *  object from which to otherwise copy construct
             */
            Base(field::FrmOnv &mbf, const Base &other) : Base(mbf) {}

            const field::FrmOnv &value() const {
                return m_mbf;
            }
        };

        class General : public Base {
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

        public:
            General(size_t nsite, size_t nelec) : Base(nsite), m_foreach(*this, nelec) {}

            General(field::FrmOnv &mbf, size_t nelec) : Base(mbf), m_foreach(*this, nelec) {}

            General(field::FrmOnv &mbf, const General &other) : General(mbf, other.m_foreach.m_nind) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
            }
        };

        class Spins : public Base {
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

            int ms2() const {
                return 2 * m_foreach.m_nind - m_mbf.m_nsite;
            }

        public:
            Spins(size_t nsite, int ms2) : Base(nsite), m_foreach(*this, ms2) {}

            Spins(field::FrmOnv &mbf, int ms2) : Base(mbf), m_foreach(*this, ms2) {}

            Spins(field::FrmOnv &mbf, const Spins &other) : Spins(mbf, other.ms2()) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
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


            size_t nelec() const {
                return m_alpha_foreach.m_nind + m_beta_foreach.m_nind;
            }

            int ms2() const {
                return m_alpha_foreach.m_nind - m_beta_foreach.m_nind;
            }

        public:
            Ms2Conserve(size_t nsite, size_t nelec, int ms2) :
                    Base(nsite),
                    m_alpha_foreach(*this, nalpha(nelec, ms2)),
                    m_beta_foreach(*this, nbeta(nelec, ms2)) {}

            Ms2Conserve(field::FrmOnv &mbf, size_t nelec, int ms2) :
                    Base(mbf),
                    m_alpha_foreach(*this, nalpha(nelec, ms2)),
                    m_beta_foreach(*this, nbeta(nelec, ms2)) {}

            Ms2Conserve(field::FrmOnv &mbf, const Ms2Conserve &other) : Ms2Conserve(mbf, other.nelec(), other.ms2()) {};

            void loop() override {
                m_alpha_foreach.loop();
            }

            size_t iiter() override {
                return m_alpha_foreach.iiter() * m_beta_foreach.m_niter + m_beta_foreach.iiter();
            }
        };
    }

    namespace bos {
        class Base : public mbf_foreach::Base {
            buffered::BosOnv m_mbf_internal;
        protected:
            field::BosOnv &m_mbf;
        public:
            Base(size_t nmode) : mbf_foreach::Base({0, nmode}), m_mbf_internal(nmode), m_mbf(m_mbf_internal) {}

            Base(field::BosOnv &mbf) : mbf_foreach::Base({0, mbf.m_nmode}), m_mbf_internal(mbf), m_mbf(mbf) {}

            Base(field::BosOnv &mbf, const Base &other) : Base(mbf) {}

            const field::BosOnv &value() const {
                return m_mbf;
            }
        };

        class General : public Base {
            struct Foreach : Unrestricted {
                General &m_context;

                Foreach(General &context, size_t nboson_max) :
                        Unrestricted(context.m_mbf.m_nmode, nboson_max), m_context(context) {}

                void body() override {
                    m_context.m_mbf = value();
                }
            };

            Foreach m_foreach;
        public:
            General(size_t nmode, size_t nboson_max) : Base(nmode), m_foreach(*this, nboson_max) {}

            General(field::BosOnv &mbf, size_t nboson_max) : Base(mbf), m_foreach(*this, nboson_max) {}

            General(field::BosOnv &mbf, const General &other) : General(mbf, other.m_foreach.m_shape[0]) {}

            void loop() override {
                m_foreach.loop();
            }

            size_t iiter() override {
                return m_foreach.iiter();
            }
        };
    }

    namespace frm_bos {

        class Base : public mbf_foreach::Base {
            buffered::FrmBosOnv m_mbf_internal;
        protected:
            field::FrmBosOnv &m_mbf;
        public:
            Base(BasisDims bd) : mbf_foreach::Base(bd), m_mbf_internal(bd), m_mbf(m_mbf_internal) {}

            Base(field::FrmBosOnv &mbf) : mbf_foreach::Base({mbf.nsite(), mbf.m_bos.m_nmode}),
                                          m_mbf_internal(mbf), m_mbf(mbf) {}

            Base(field::FrmBosOnv &mbf, const Base &other) : Base(mbf) {}

            const field::FrmBosOnv &value() const {
                return m_mbf;
            }
        };

        template<typename frm_foreach_t, typename bos_foreach_t>
        class Product : public Base {
            static_assert(std::is_base_of<mbf_foreach::frm::Base, frm_foreach_t>::value,
                          "template arg must be derived from the base type of fermion foreach iterators");
            static_assert(std::is_base_of<mbf_foreach::bos::Base, bos_foreach_t>::value,
                          "template arg must be derived from the base type of boson foreach iterators");

            struct BosForeach : public bos_foreach_t {
                BosForeach(field::BosOnv &mbf, const bos_foreach_t &other) : bos_foreach_t(mbf, other) {}
                void body() override {

                }
            };

            struct FrmForeach : public frm_foreach_t {
                FrmForeach(field::FrmOnv &mbf, const frm_foreach_t &other) : frm_foreach_t(mbf, other) {}
                void body() override {

                }
            };

        public:
            FrmForeach m_frm_foreach;
            BosForeach m_bos_foreach;

            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) :
                    Base({static_cast<const mbf_foreach::frm::Base &>(frm_foreach).m_bd.m_nsite,
                          static_cast<const mbf_foreach::bos::Base &>(bos_foreach).m_bd.m_nmode}),
                    m_frm_foreach(m_mbf.m_frm, frm_foreach), m_bos_foreach(m_mbf.m_bos, bos_foreach) {}

//            Product(field::FrmBosOnv &mbf, const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) :
//                    Base(mbf), m_frm_foreach(mbf.m_frm, frm_foreach), m_bos_foreach(mbf.m_bos, bos_foreach){}
//
//            Product(field::FrmBosOnv &mbf, const Product &other) :
//                    Product(mbf, other.m_frm_foreach, other.m_bos_foreach){}

            void body() override {

            }

            void loop() override {

            }

            size_t iiter() override {
                return 0;
            }
        };
    }

}

#endif //M7_MBFFOREACH_H
