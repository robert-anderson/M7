//
// Created by rja on 27/03/2022.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H

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

        Base(BasisData bd, size_t niter) : m_bd(bd), m_niter(niter) {}

        virtual ~Base() {}

        template<typename mbf_t>
        using function_t = std::function<void(const mbf_t &)>;
    protected:
        virtual void frm_loop(const function_t<field::FrmOnv> &fn) {};

        virtual void bos_loop(const function_t<field::BosOnv> &fn) {};

        virtual void frmbos_loop(const function_t<field::FrmBosOnv> &fn) {};


    public:
        void loop(const function_t<field::FrmOnv> &fn) { frm_loop(fn); }

        void loop(const function_t<field::BosOnv> &fn) { bos_loop(fn); }

        void loop(const function_t<field::FrmBosOnv> &fn) { frmbos_loop(fn); }
    };

    struct PairBase {

        template<typename mbf_t>
        using function_t = std::function<void(const mbf_t &, size_t, const mbf_t &, size_t)>;
    protected:
        virtual void frm_loop(const function_t<field::FrmOnv> &fn) {};

        virtual void bos_loop(const function_t<field::BosOnv> &fn) {};

        virtual void frmbos_loop(const function_t<field::FrmBosOnv> &fn) {};


    public:
        void loop(const function_t<field::FrmOnv> &fn) { frm_loop(fn); }

        void loop(const function_t<field::BosOnv> &fn) { bos_loop(fn); }

        void loop(const function_t<field::FrmBosOnv> &fn) { frmbos_loop(fn); }
    };

    template<typename mbf_t, typename foreach_t>
    struct Pair : PairBase {
        static_assert(std::is_base_of<Base, foreach_t>::value,
                      "template arg must be derived from mbf_foreach::Base");
        foreach_t m_foreach_outer;
        foreach_t m_foreach_inner;

        Pair(const foreach_t &foreach) : m_foreach_outer(foreach), m_foreach_inner(foreach) {}

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

            Base(size_t nsite, size_t niter) : mbf_foreach::Base({nsite, 0ul}, niter), m_mbf(nsite) {}
        };

        struct NumberConserve : Base {
            const size_t m_nelec;

            NumberConserve(size_t nsite, size_t nelec, size_t niter) : Base(nsite, niter), m_nelec(nelec) {}
        };

        struct General : NumberConserve {
            basic_foreach::rtnd::Ordered<true, true> m_foreach;

            General(size_t nsite, size_t nelec) : NumberConserve(nsite, nelec, 1),
                                                  m_foreach(m_bd.m_nspinorb, m_nelec) {}

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
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override { loop(fn); }
        };


        struct Spins : NumberConserve {
            const int m_ms2;
            basic_foreach::rtnd::Ordered<true, true> m_foreach;

            Spins(size_t nsite, int ms2) : NumberConserve(nsite, nsite + ms2 / 2, 1), m_ms2(ms2),
                                           m_foreach(m_bd.m_nsite, (m_bd.m_nsite + m_ms2) / 2) {}

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
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override { loop(fn); }
        };


        struct Ms2Conserve : NumberConserve {
            const int m_ms2;
            basic_foreach::rtnd::Ordered<true, true> m_alpha_foreach;
            basic_foreach::rtnd::Ordered<true, true> m_beta_foreach;

            Ms2Conserve(size_t nsite, size_t nelec, int ms2) : NumberConserve(nsite, nelec, 1), m_ms2(ms2),
                                                               m_alpha_foreach(m_bd.m_nsite,
                                                                               ci_utils::nalpha(m_nelec, m_ms2)),
                                                               m_beta_foreach(m_bd.m_nsite,
                                                                              ci_utils::nbeta(m_nelec, m_ms2)) {}

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
            void frm_loop(const std::function<void(const field::FrmOnv &)> &fn) override { loop(fn); }

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

            Base(size_t nmode, size_t niter) : mbf_foreach::Base({0ul, nmode}, niter), m_mbf(nmode) {}
        };

        struct NumberConserve : Base {
            const size_t m_nboson;

            NumberConserve(size_t nmode, size_t nboson, size_t niter) : Base(nmode, niter), m_nboson(nboson) {}
        };

        struct GeneralClosed : NumberConserve {
            basic_foreach::rtnd::Ordered<false, true> m_foreach;

            GeneralClosed(size_t nmode, size_t nboson) : NumberConserve(nmode, nboson, 1),
                                                         m_foreach(nmode, nboson) {}

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
            void bos_loop(const std::function<void(const field::BosOnv &)> &fn) override { loop(fn); }

        };

        struct GeneralOpen : Base {
            basic_foreach::rtnd::Unrestricted m_foreach;

            GeneralOpen(size_t nmode, size_t nboson_max) : Base(nmode, 1),
                                                           m_foreach(nmode, nboson_max + 1) {}

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
            void bos_loop(const std::function<void(const field::BosOnv &)> &fn) override { loop(fn); }

        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::BosOnv, foreach_t> {
            static_assert(std::is_base_of<bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::bos::Base");
            typedef mbf_foreach::Pair<field::BosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void bos_loop(const mbf_foreach::PairBase::function_t<field::BosOnv> &fn) override {base_t::loop(fn);}
        };
    }

    namespace frm_bos {
        struct Base : mbf_foreach::Base {
            buffered::FrmBosOnv m_mbf;
            Base(BasisData bd, size_t niter) : mbf_foreach::Base(bd, niter), m_mbf(bd) {}
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
            static std::pair<const frm::Base&, const bos::Base&> cast_to_bases(
                    const frm_foreach_t& frm_foreach, const bos_foreach_t& bos_foreach) {
                return {frm_foreach, bos_foreach};
            }
            static BasisData make_bd(const frm_foreach_t& frm_foreach, const bos_foreach_t& bos_foreach) {
                auto bases = cast_to_bases(frm_foreach, bos_foreach);
                return {bases.first.m_bd.m_nsite, bases.second.m_bd.m_nmode};
            }
            static size_t make_niter(const frm_foreach_t& frm_foreach, const bos_foreach_t& bos_foreach) {
                auto bases = cast_to_bases(frm_foreach, bos_foreach);
                return bases.first.m_niter * bases.second.m_niter;
            }

        public:
            Product(const frm_foreach_t& frm_foreach, const bos_foreach_t& bos_foreach):
                    Base(make_bd(frm_foreach, bos_foreach), make_niter(frm_foreach, bos_foreach)),
                    m_frm_foreach(frm_foreach), m_bos_foreach(bos_foreach){}


            template<typename fn_t>
            void loop(const fn_t &fn) {
                functor_utils::assert_prototype<void(const field::FrmBosOnv &)>(fn);

                auto frm_loop_fn = [this, &fn](const field::FrmOnv& frm_mbf){
                    m_mbf.m_frm = frm_mbf;
                    auto bos_loop_fn = [this, &fn](const field::BosOnv& bos_mbf) {
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


        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::FrmBosOnv, foreach_t> {
            static_assert(std::is_base_of<frm_bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm_bos::Base");
            typedef mbf_foreach::Pair<field::FrmBosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void frmbos_loop(const mbf_foreach::PairBase::function_t<field::FrmBosOnv> &fn) override {base_t::loop(fn);}
        };
    }

}


#endif //M7_MBFFOREACH_H
