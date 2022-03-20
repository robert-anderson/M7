//
// Created by anderson on 2/9/22.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H


#include <utility>

#include <M7_lib/field/Fields.h>
#include <M7_lib/foreach/ForeachVirtual.h>
#include <M7_lib/table/BufferedFields.h>

namespace mbf_foreach {
    /**
     * all MBFs have runtime-determined dimensions
     */
    using namespace foreach_virtual::rtnd;

    /**
     * untemplated, polymorphic base class
     */
    struct MbfForeach {
        const BasisData m_bd;

        MbfForeach(BasisData bd) : m_bd(bd) {}
        MbfForeach(const MbfForeach &other) : m_bd(other.m_bd) {}
        virtual ~MbfForeach(){}

        /**
         * iterates over all values and iiters
         */
        void loop() {
            try {throwing_loop();}
            catch (const ExitLoop&){}
        }

        virtual void throwing_loop() = 0;
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
        Base(BasisData bd, body_fn_t body_fn = {}, field_t *mbf = nullptr) :
                MbfForeach(bd), m_mbf_internal(bd), m_mbf(mbf ? mbf : &m_mbf_internal), m_body_fn(std::move(body_fn)) {}

        Base(const Base &other, field_t *mbf = nullptr) : Base(other.m_bd, other.m_body_fn, mbf) {}

        virtual void body() {
            if (m_body_fn) m_body_fn(*m_mbf, iiter());
        }

        const field_t& mbf() const {
            return *m_mbf;
        }
    };

    namespace frm {

        class Base : public mbf_foreach::Base<defs::Frm> {
        public:
            static constexpr size_t mbf_ind = defs::Frm;
            Base(size_t nsite, body_fn_t body_fn = {}, field_t *mbf = nullptr);

            Base(const mbf_foreach::frm::Base &other, field_t *mbf);
        };

        class General : public Base {
            struct Foreach : Ordered<> {
                General &m_context;

                Foreach(General &context, size_t nelec);

                void body(const inds_t &value, size_t iiter) override;
            };

            Foreach m_foreach;

        public:
            General(size_t nsite, size_t nelec, body_fn_t body_fn = {}, field_t *mbf = nullptr);

            General(const General &other, field_t *mbf = nullptr);

            void throwing_loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };

        class Spins : public Base {
            struct Foreach : Ordered<> {
                Spins &m_context;

                Foreach(Spins &context, int ms2);

                void body(const inds_t &value, size_t iiter) override;

            };

            Foreach m_foreach;

            int ms2() const;

        public:
            Spins(size_t nsite, int ms2, body_fn_t body_fn = {}, field_t *mbf = nullptr);

            Spins(const Spins &other, field_t *mbf = nullptr);

            void throwing_loop() override;

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

                void body(const inds_t &value, size_t iiter) override;
            };

            /**
             * iterator over outer spin channel
             */
            struct AlphaForeach : Foreach {
                AlphaForeach(Ms2Conserve &context, size_t nelec) : Foreach(context, nelec) {}

                void body(const inds_t &value, size_t iiter) override;
            };

            AlphaForeach m_alpha_foreach;
            BetaForeach m_beta_foreach;


            size_t nelec() const;

            int ms2() const;

        public:
            Ms2Conserve(size_t nsite, size_t nelec, int ms2, body_fn_t body_fn = {}, field_t *mbf = nullptr);

            Ms2Conserve(const Ms2Conserve &other, field_t *mbf = nullptr);

            void throwing_loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };
    }

    namespace bos {

        class Base : public mbf_foreach::Base<defs::Bos> {
        public:
            static constexpr size_t mbf_ind = defs::Bos;
            Base(size_t nmode, body_fn_t body_fn = {}, field_t *mbf = nullptr);

            Base(const Base &other, field_t *mbf = nullptr);
        };

        class GeneralOpen : public Base {
            struct Foreach : Unrestricted {
                GeneralOpen &m_context;

                Foreach(GeneralOpen &context, size_t nboson_max);

                void body(const inds_t &value, size_t iiter) override;

                size_t nboson_max() const;
            };

            Foreach m_foreach;
        public:
            GeneralOpen(size_t nmode, size_t nboson_max, body_fn_t body_fn = {}, field::BosOnv *mbf = nullptr);

            GeneralOpen(const GeneralOpen &other, field::BosOnv *mbf = nullptr);

            void throwing_loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };

        class GeneralClosed : public Base {
            struct Foreach : Ordered<false, true> {
                GeneralClosed &m_context;

                Foreach(GeneralClosed &context, size_t nboson):
                    Ordered<false, true>(context.m_bd.m_nmode, nboson), m_context(context){}

                void body(const inds_t &value, size_t iiter) override {
                    m_context.m_mbf->set_ops(value);
                    m_context.body();
                }
            };

            Foreach m_foreach;

        public:
            GeneralClosed(size_t nmode, size_t nboson, body_fn_t body_fn = {}, field::BosOnv *mbf = nullptr):
                    Base(nmode, body_fn, mbf), m_foreach(*this, nboson) {}

            GeneralClosed(const GeneralClosed &other, field::BosOnv *mbf = nullptr):
                    GeneralClosed(other.m_bd.m_nmode, other.m_foreach.m_nind, other.m_body_fn, mbf) {}

            void throwing_loop() override;

            size_t iiter() const override;

            size_t niter() const override;
        };
    }

    namespace frm_bos {

        class Base : public mbf_foreach::Base<defs::FrmBos> {
        public:
            static constexpr size_t mbf_ind = defs::FrmBos;
            Base(BasisData bd, body_fn_t body_fn = {}, field_t *mbf = nullptr) :
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
                    m_context.m_bos_foreach.throwing_loop();
                }
            };

        public:
            FrmForeach m_frm_foreach;
            BosForeach m_bos_foreach;

            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach,
                        body_fn_t body_fn = {}, field_t *mbf=nullptr):
                    Base({frm_foreach.m_bd.m_nsite, bos_foreach.m_bd.m_nmode}, std::move(body_fn), mbf),
                    m_frm_foreach(*this, frm_foreach), m_bos_foreach(*this, bos_foreach) {}

            Product(const Product &other, field_t *mbf = nullptr) :
                    Product(other.m_frm_foreach, other.m_bos_foreach, other.m_body_fn, mbf) {}

            void throwing_loop() override {
                m_frm_foreach.throwing_loop();
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

        /**
         * convenient definition for case when the boson sector iteration is "general" i.e. unconstrained and boson
         * number non-conserving upto a maximum occupation cutoff
         */
        template<typename frm_foreach_t>
        using BosGeneralOpen = Product<frm_foreach_t, bos::GeneralOpen>;
    }

    class PairBase : public MbfForeach {
    public:
        PairBase(BasisData bd): MbfForeach(bd){}
        PairBase(const MbfForeach& other): MbfForeach(other){}
        virtual ~PairBase(){}
        virtual size_t nrow() = 0;
    };

    template<typename foreach_t>
    class Pair : public PairBase {
        static constexpr size_t mbf_ind = foreach_t::mbf_ind;
        static_assert(std::is_base_of<Base<mbf_ind>, foreach_t>::value,
                "template arg foreach_t is not compatible with mbf_ind");
    protected:
        typedef field::mbf_t<mbf_ind> field_t;
        typedef std::function<void(const field_t&, size_t, const field_t&, size_t)> body_fn_t;

        body_fn_t m_body_fn;

        struct Inner : public foreach_t {
            Pair& m_context;
            Inner(Pair& context, const foreach_t &other) : foreach_t(other, nullptr), m_context(context){}
            void body() override {
                foreach_t::body();
                m_context.body();
            }
        };

        struct Outer : public foreach_t {
            Pair& m_context;
            Outer(Pair& context, const foreach_t &other) : foreach_t(other, nullptr), m_context(context){}
            void body() override {
                foreach_t::body();
                m_context.m_inner.throwing_loop();
            }
        };

        Inner m_inner;
        Outer m_outer;

        void body() {
            if (m_body_fn) m_body_fn(m_outer.mbf(), m_outer.iiter(), m_inner.mbf(), m_inner.iiter());
        }
    public:
        void throwing_loop() override {
            m_outer.throwing_loop();
        }

        size_t iiter() const override {
            return m_outer.iiter()*m_inner.niter()+m_inner.iiter();
        }

        size_t niter() const override {
            return m_outer.niter()*m_inner.niter();
        }

    private:
        size_t nrow() override {
            return m_outer.niter();
        }

    public:

        Pair(const foreach_t &foreach, body_fn_t body_fn = {}):
                PairBase(foreach), m_body_fn(body_fn), m_inner(*this, foreach), m_outer(*this, foreach){}

    };
}

#endif //M7_MBFFOREACH_H
