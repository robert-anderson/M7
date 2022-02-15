//
// Created by anderson on 2/9/22.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H


#include <src/core/field/Fields.h>
#include <src/core/util/ForeachVirtual.h>
#include <src/core/table/BufferedFields.h>
#include <utility>

namespace mbf_foreach {
    /**
     * all MBFs have runtime-determined dimensions
     */
    using namespace foreach_virtual::rtnd;

    /**
     * untemplated, polymorphic base class
     */
    struct MbfForeach {
        const BasisDims m_bd;

        MbfForeach(BasisDims bd) : m_bd(bd) {}

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

    template<size_t mbf_ind>
    class Base : public MbfForeach {
        typedef buffered::mbf_t<mbf_ind> buffered_t;
        buffered_t m_mbf_internal;
    protected:
        typedef field::mbf_t<mbf_ind> field_t;
        field_t *m_mbf;

        typedef std::function<void(const field_t&, size_t)> body_fn_t;
        body_fn_t m_body_fn;
    public:
        Base(BasisDims bd, body_fn_t body_fn, field_t *mbf = nullptr) :
                MbfForeach(bd), m_mbf_internal(bd), m_mbf(mbf ? mbf : &m_mbf_internal), m_body_fn(std::move(body_fn)) {}

        Base(const Base &other, field_t *mbf = nullptr) : Base(other.m_bd, other.m_body_fn, mbf) {}

        virtual void body() {
            if (m_body_fn) m_body_fn(*m_mbf, iiter());
        }
    };

    namespace frm {

        class Base : public mbf_foreach::Base<defs::Frm> {
        public:
            Base(size_t nsite, body_fn_t body_fn, field_t *mbf = nullptr);

            Base(const mbf_foreach::frm::Base &other, field_t *mbf);
        };

        class General : public Base {
            struct Foreach : Ordered<> {
                General &m_context;

                Foreach(General &context, size_t nelec);

                void body() override;
            };

            Foreach m_foreach;

        public:
            General(size_t nsite, size_t nelec, body_fn_t body_fn, field_t *mbf = nullptr);

            General(const General &other, field_t *mbf = nullptr);

            void loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };

        class Spins : public Base {
            struct Foreach : Ordered<> {
                Spins &m_context;

                Foreach(Spins &context, int ms2);

                void body() override;

            };

            Foreach m_foreach;

            int ms2() const;

        public:
            Spins(size_t nsite, int ms2, body_fn_t body_fn, field_t *mbf = nullptr);

            Spins(const Spins &other, field_t *mbf = nullptr);

            void loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };

        using namespace ci_utils;
        class Ms2Conserve : public Base {
            /**
             * iterator over a single spin channel
             */
            struct Foreach : Ordered<> {
                Ms2Conserve &m_context;

                Foreach(Ms2Conserve &context, size_t nelec);
            };

            /**
             * iterator over inner spin channel
             */
            struct BetaForeach : Foreach {
                BetaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body() override;
            };

            /**
             * iterator over outer spin channel
             */
            struct AlphaForeach : Foreach {
                AlphaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body() override;
            };

            AlphaForeach m_alpha_foreach;
            BetaForeach m_beta_foreach;


            size_t nelec() const;

            int ms2() const;

        public:
            Ms2Conserve(size_t nsite, size_t nelec, int ms2, body_fn_t body_fn, field_t *mbf = nullptr);

            Ms2Conserve(const Ms2Conserve &other, field_t *mbf = nullptr);

            void loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };
    }

    namespace bos {

        class Base : public mbf_foreach::Base<defs::Bos> {
        public:
            Base(size_t nmode, body_fn_t body_fn, field_t *mbf = nullptr);

            Base(const Base &other, field_t *mbf = nullptr);
        };

        class General : public Base {
            struct Foreach : Unrestricted {
                General &m_context;

                Foreach(General &context, size_t nboson_max);

                void body() override;

                size_t nboson_max() const;
            };

            Foreach m_foreach;
        public:
            General(size_t nmode, size_t nboson_max, body_fn_t body_fn, field::BosOnv *mbf = nullptr);

            General(const General &other, field::BosOnv *mbf = nullptr);

            void loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };
    }

    namespace frm_bos {

        class Base : public mbf_foreach::Base<defs::FrmBos> {
        public:
            Base(BasisDims bd, body_fn_t body_fn, field_t *mbf = nullptr) :
                    mbf_foreach::Base<defs::FrmBos>(bd, std::move(body_fn), mbf) {}

            Base(const Base &other, field_t *mbf = nullptr) :
                    Base(other.m_bd, other.m_body_fn, mbf) {}
        };

        template<typename frm_foreach_t, typename bos_foreach_t>
        class Product : public Base {
            static_assert(std::is_base_of<frm::Base, frm_foreach_t>::value,
                          "template arg must be derived from the base type of fermion foreach iterators");
            static_assert(std::is_base_of<bos::Base, bos_foreach_t>::value,
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

            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach,
                        body_fn_t body_fn, field_t *mbf=nullptr):
                    Base({frm_foreach.m_bd.m_nsite, bos_foreach.m_bd.m_nmode}, std::move(body_fn), mbf),
                    m_frm_foreach(*this, frm_foreach), m_bos_foreach(*this, bos_foreach) {}

            Product(const Product &other, field_t *mbf = nullptr) :
                    Base(other.m_frm_foreach, other.m_bos_foreach, other.m_body_fn, mbf) {}

            void loop() override {
                m_frm_foreach.loop();
            }

            size_t iiter() const override {
                auto iiter_frm = m_frm_foreach.iiter();
                auto iiter_bos = m_bos_foreach.iiter();
                auto niter_bos = m_bos_foreach.niter();
                return iiter_frm * niter_bos + iiter_bos;
            }

            size_t niter() const override {
                auto niter_frm = m_frm_foreach.niter();
                auto niter_bos = m_bos_foreach.niter();
                return niter_frm * niter_bos;
            }
        };
    }
}

#endif //M7_MBFFOREACH_H
