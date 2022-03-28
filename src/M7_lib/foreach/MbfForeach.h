//
// Created by rja on 27/03/2022.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H

#include <utility>

#include "BasicForeach.h"
#include "M7_lib/field/Fields.h"
#include "M7_lib/table/BufferedFields.h"

namespace mbf_foreach {
    /**
     * all MBFs have runtime-determined dimensions
     */
    using namespace basic_foreach::rtnd;

    /**
     * untemplated, polymorphic base class
     */
    struct Base {
        const BasisData m_bd;
        const size_t m_niter;

        Base(BasisData bd, size_t niter);

        virtual ~Base() {}

        template<typename mbf_t>
        using function_t = std::function<void(const mbf_t &)>;
    protected:
        virtual void frm_loop(const function_t<field::FrmOnv> &fn) {};

        virtual void bos_loop(const function_t<field::BosOnv> &fn) {};

        virtual void frmbos_loop(const function_t<field::FrmBosOnv> &fn) {};


    public:
        void loop(const function_t<field::FrmOnv> &fn);

        void loop(const function_t<field::BosOnv> &fn);

        void loop(const function_t<field::FrmBosOnv> &fn);
    };

    struct PairBase {

        const size_t m_nrow;
        const size_t m_niter;

        template<typename mbf_t>
        using function_t = std::function<void(const mbf_t &, size_t, const mbf_t &, size_t)>;
    protected:
        virtual void frm_loop(const function_t<field::FrmOnv> &fn) {};

        virtual void bos_loop(const function_t<field::BosOnv> &fn) {};

        virtual void frmbos_loop(const function_t<field::FrmBosOnv> &fn) {};


    public:
        PairBase(size_t nrow);

        virtual ~PairBase() {}

        void loop(const function_t<field::FrmOnv> &fn);

        void loop(const function_t<field::BosOnv> &fn);

        void loop(const function_t<field::FrmBosOnv> &fn);
    };

    template<typename mbf_t, typename foreach_t>
    struct Pair : PairBase {
        static_assert(std::is_base_of<Base, foreach_t>::value,
                      "template arg must be derived from mbf_foreach::Base");
        foreach_t m_foreach_outer;
        foreach_t m_foreach_inner;

        Pair(const foreach_t &foreach) :
                PairBase(static_cast<const mbf_foreach::Base &>(foreach).m_niter),
                m_foreach_outer(foreach), m_foreach_inner(foreach) {}

        template<typename fn_t>
        void loop(const fn_t &fn) {
            functor_utils::assert_prototype<void(const mbf_t &, size_t, const mbf_t &, size_t)>(fn);
            size_t iouter = 0ul;
            auto outer_fn = [this, &fn, &iouter](const mbf_t &outer) {
                size_t iinner = 0ul;
                auto inner_fn = [this, &fn, &outer, &iouter, &iinner](const mbf_t &inner) {
                    fn(outer, iouter, inner, iinner);
                    ++iinner;
                };
                m_foreach_inner.loop(inner_fn);
                ++iouter;
            };
            m_foreach_outer.loop(outer_fn);
        }
    };

    namespace frm {

        struct Base : mbf_foreach::Base {
            buffered::FrmOnv m_mbf;

            Base(size_t nsite, size_t niter);
        };

        struct NumberConserve : Base {
            const size_t m_nelec;

            NumberConserve(size_t nsite, size_t nelec, size_t niter);
        };

        struct General : NumberConserve {
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_foreach;

            General(size_t nsite, size_t nelec);

        public:
            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::FrmOnv &)>(fn);
                auto loop_fn = [this, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    m_mbf = inds;
                    fn(m_mbf);
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override;
        };


        struct Spins : NumberConserve {
            const int m_ms2;
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_foreach;

            Spins(size_t nsite, int ms2);

        public:
            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::FrmOnv &)>(fn);
                auto loop_fn = [this, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    m_mbf.set_spins(inds);
                    fn(m_mbf);
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override;
        };


        struct Ms2Conserve : NumberConserve {
            const int m_ms2;
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_alpha_foreach, m_beta_foreach;

            static size_t niter(size_t nsite, size_t nelec, int ms2);

            Ms2Conserve(size_t nsite, size_t nelec, int ms2);

        public:
            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::FrmOnv &)>(fn);
                auto alpha_fn = [this, &fn](const basic_foreach::rtnd::inds_t &alpha_inds) {
                    m_mbf.put_spin_channel(0, false);
                    m_mbf.set(0, alpha_inds);
                    auto beta_fn = [this, &fn](const basic_foreach::rtnd::inds_t &beta_inds) {
                        m_mbf.put_spin_channel(1, false);
                        m_mbf.set(m_mbf.m_nsite, beta_inds);
                        fn(m_mbf);
                    };
                    m_beta_foreach.loop(beta_fn);
                };
                m_alpha_foreach.loop(alpha_fn);
            }

        protected:
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override;

        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::FrmOnv, foreach_t> {
            static_assert(std::is_base_of<frm::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm::Base");
            typedef mbf_foreach::Pair<field::FrmOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void frm_loop(const mbf_foreach::PairBase::function_t<field::FrmOnv> &fn) override {
                base_t::loop(fn);
            }
        };
    }

    namespace bos {
        struct Base : mbf_foreach::Base {
            buffered::BosOnv m_mbf;

            Base(size_t nmode, size_t niter);
        };

        struct NumberConserve : Base {
            const size_t m_nboson;

            NumberConserve(size_t nmode, size_t nboson, size_t niter);
        };

        struct GeneralClosed : NumberConserve {
            typedef basic_foreach::rtnd::Ordered<false, true> foreach_t;
            foreach_t m_foreach;

            GeneralClosed(size_t nmode, size_t nboson);

        public:
            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::BosOnv &)>(fn);
                auto loop_fn = [this, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    m_mbf.set_ops(inds);
                    fn(m_mbf);
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void bos_loop(const std::function<void(const field::BosOnv &)> &fn) override;

        };

        struct GeneralOpen : Base {
            typedef basic_foreach::rtnd::Unrestricted foreach_t;
            foreach_t m_foreach;

            GeneralOpen(size_t nmode, size_t nboson_max);

        public:
            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::BosOnv &)>(fn);
                auto loop_fn = [this, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    m_mbf = inds;
                    fn(m_mbf);
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void bos_loop(const std::function<void(const field::BosOnv &)> &fn) override;

        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::BosOnv, foreach_t> {
            static_assert(std::is_base_of<bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::bos::Base");
            typedef mbf_foreach::Pair<field::BosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void bos_loop(const mbf_foreach::PairBase::function_t<field::BosOnv> &fn) override;
        };

        template<typename foreach_t>
        void Pair<foreach_t>::bos_loop(const PairBase::function_t<field::BosOnv> &fn) { base_t::loop(fn); }
    }

    namespace frm_bos {
        struct Base : mbf_foreach::Base {
            buffered::FrmBosOnv m_mbf;

            Base(BasisData bd, size_t niter);
        };

        template<typename frm_foreach_t, typename bos_foreach_t>
        struct Product : Base {
            static_assert(std::is_base_of<frm::Base, frm_foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm::Base");
            static_assert(std::is_base_of<bos::Base, bos_foreach_t>::value,
                          "template arg must be derived from mbf_foreach::bos::Base");
            frm_foreach_t m_frm_foreach;
            bos_foreach_t m_bos_foreach;

        private:
            static std::pair<const frm::Base &, const bos::Base &> cast_to_bases(
                    const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) {
                return {frm_foreach, bos_foreach};
            }

            static BasisData make_bd(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) {
                auto bases = cast_to_bases(frm_foreach, bos_foreach);
                return {bases.first.m_bd.m_nsite, bases.second.m_bd.m_nmode};
            }

            static size_t make_niter(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) {
                auto bases = cast_to_bases(frm_foreach, bos_foreach);
                return bases.first.m_niter * bases.second.m_niter;
            }

        public:
            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) :
                    Base(make_bd(frm_foreach, bos_foreach), make_niter(frm_foreach, bos_foreach)),
                    m_frm_foreach(frm_foreach), m_bos_foreach(bos_foreach) {}


            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::FrmBosOnv &)>(fn);

                auto frm_loop_fn = [this, &fn](const field::FrmOnv &frm_mbf) {
                    m_mbf.m_frm = frm_mbf;
                    auto bos_loop_fn = [this, &fn](const field::BosOnv &bos_mbf) {
                        m_mbf.m_bos = bos_mbf;
                        fn(m_mbf);
                    };
                    m_bos_foreach.loop(bos_loop_fn);
                };
                m_frm_foreach.loop(frm_loop_fn);
            }

        protected:
            void frmbos_loop(const std::function<void(const field::FrmBosOnv &)> &fn) override { loop(fn); }
        };

        /**
         * convenient partial specialization for fermions coupled to a closed bosonic quantum system
         * @tparam frm_foreach_t
         *  fermion foreach iterator type
         */
        template<typename frm_foreach_t>
        struct ClosedProduct : Product<frm_foreach_t, bos::GeneralClosed> {
            ClosedProduct(const frm_foreach_t &frm_foreach, size_t nmode, size_t nboson) :
                    Product<frm_foreach_t, bos::GeneralClosed>(frm_foreach, {nmode, nboson}) {}
        };

        /**
         * convenient partial specialization for fermions coupled to an open bosonic quantum system
         * @tparam frm_foreach_t
         *  fermion foreach iterator type
         */
        template<typename frm_foreach_t>
        struct OpenProduct : Product<frm_foreach_t, bos::GeneralOpen> {
            OpenProduct(const frm_foreach_t &frm_foreach, size_t nmode, size_t nboson_max) :
                    Product<frm_foreach_t, bos::GeneralOpen>(frm_foreach, {nmode, nboson_max}) {}
        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::FrmBosOnv, foreach_t> {
            static_assert(std::is_base_of<frm_bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm_bos::Base");
            typedef mbf_foreach::Pair<field::FrmBosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void frmbos_loop(const mbf_foreach::PairBase::function_t<field::FrmBosOnv> &fn) override {
                base_t::loop(fn);
            }
        };

    }

}


#endif //M7_MBFFOREACH_H
