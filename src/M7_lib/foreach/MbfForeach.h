//
// Created by Robert J. Anderson on 27/03/2022.
//

#ifndef M7_MBFFOREACH_H
#define M7_MBFFOREACH_H

#include <utility>

#include "BasicForeach.h"
#include "M7_lib/field/Fields.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/basis/Suites.h"

namespace mbf_foreach {
    /**
     * all MBFs have runtime-determined dimensions
     */
    using namespace basic_foreach::rtnd;

    /**
     * untemplated, polymorphic base class
     */
    struct Base {
        using prototype_t = void();
        using function_t = std::function<prototype_t>;

        const size_t m_niter;

        Base(size_t niter);

        virtual ~Base() {}

    protected:
        virtual void frm_loop(field::FrmOnv &mbf, const function_t &fn) {};

        virtual void bos_loop(field::BosOnv &mbf, const function_t &fn) {};

        virtual void frmbos_loop(field::FrmBosOnv &mbf, const function_t &fn) {};


    public:
        template<typename fn_t>
        void loop(field::FrmOnv &mbf, const fn_t &fn) {
            functor_utils::assert_prototype<void()>(fn);
            frm_loop(mbf, fn);
        }

        template<typename fn_t>
        void loop(field::BosOnv &mbf, const fn_t &fn) {
            functor_utils::assert_prototype<void()>(fn);
            bos_loop(mbf, fn);
        }

        template<typename fn_t>
        void loop(field::FrmBosOnv &mbf, const fn_t &fn) {
            functor_utils::assert_prototype<void()>(fn);
            frmbos_loop(mbf, fn);
        }
    };

    struct PairBase {
        using prototype_t = void(size_t, size_t);
        using function_t = std::function<prototype_t>;

        const size_t m_nrow;
        const size_t m_niter;

    protected:
        virtual void frm_loop(field::FrmOnv &bra, field::FrmOnv &ket, const function_t &fn) {};
        virtual void bos_loop(field::BosOnv &bra, field::BosOnv &ket, const function_t &fn) {};
        virtual void frmbos_loop(field::FrmBosOnv &bra, field::FrmBosOnv &ket, const function_t &fn) {};

    public:
        PairBase(size_t nrow);

        virtual ~PairBase() {}

        template<typename fn_t>
        void loop(field::FrmOnv &bra, field::FrmOnv &ket, const fn_t &fn) {
            functor_utils::assert_prototype<prototype_t>(fn);
            frm_loop(bra, ket, fn);
        }

        template<typename fn_t>
        void loop(field::BosOnv &bra, field::BosOnv &ket, const fn_t &fn) {
            functor_utils::assert_prototype<prototype_t>(fn);
            bos_loop(bra, ket, fn);
        }

        template<typename fn_t>
        void loop(field::FrmBosOnv &bra, field::FrmBosOnv &ket, const fn_t &fn) {
            functor_utils::assert_prototype<prototype_t>(fn);
            frmbos_loop(bra, ket, fn);
        }
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
        void loop_fn(mbf_t& outer, mbf_t& inner, const fn_t &fn) {
            functor_utils::assert_prototype<prototype_t>(fn);
            size_t iouter = 0ul;
            auto outer_fn = [this, &outer, &inner, &fn, &iouter]() {
                size_t iinner = 0ul;
                auto inner_fn = [&fn, &iouter, &iinner]() {
                    fn(iouter, iinner);
                    ++iinner;
                };
                m_foreach_inner.loop_fn(inner, inner_fn);
                ++iouter;
            };
            m_foreach_outer.loop_fn(outer, outer_fn);
        }
    };

    namespace frm {

        struct Base : mbf_foreach::Base {
            const sys::frm::Sector m_sector;
            Base(const sys::frm::Sector& sector, size_t niter);
        protected:
            void verify_mbf(const field::FrmOnv& mbf) {
                REQUIRE_TRUE(mbf.m_basis==m_sector.m_basis, "given MBF has incorrect basis");
            }
        };

        struct General : Base {
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_foreach;

        public:

            General(const sys::frm::Sector& sector);

            template<typename fn_t>
            void loop_fn(field::FrmOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                verify_mbf(mbf);
                auto loop_fn = [&mbf, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    mbf = inds;
                    fn();
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void frm_loop(field::FrmOnv &mbf, const function_t &fn) override;
        };


        struct Spins : Base {
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_foreach;

            Spins(const sys::frm::Sector& sector);

        public:
            template<typename fn_t>
            void loop_fn(field::FrmOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                verify_mbf(mbf);
                auto loop_fn = [&mbf, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    mbf.set_spins(inds);
                    fn();
                };
                m_foreach.loop(loop_fn);
            }

        protected:
            void frm_loop(field::FrmOnv &mbf, const function_t &fn) override;
        };


        struct Ms2Conserve : Base {
            typedef basic_foreach::rtnd::Ordered<true, true> foreach_t;
            foreach_t m_alpha_foreach, m_beta_foreach;

            static size_t niter(const sys::frm::Sector& sector);

            Ms2Conserve(const sys::frm::Sector& sector);

        public:
            template<typename fn_t>
            void loop_fn(field::FrmOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                verify_mbf(mbf);
                auto alpha_fn = [this, &mbf, &fn](const basic_foreach::rtnd::inds_t &alpha_inds) {
                    mbf.put_spin_channel(0, false);
                    mbf.set(0, alpha_inds);
                    auto beta_fn = [&mbf, &fn](const basic_foreach::rtnd::inds_t &beta_inds) {
                        mbf.put_spin_channel(1, false);
                        mbf.set(mbf.m_basis.m_nsite, beta_inds);
                        fn();
                    };
                    m_beta_foreach.loop(beta_fn);
                };
                m_alpha_foreach.loop(alpha_fn);
            }

        protected:
            void frm_loop(field::FrmOnv &mbf, const function_t &fn) override;

        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::FrmOnv, foreach_t> {
            static_assert(std::is_base_of<frm::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm::Base");
            typedef mbf_foreach::Pair<field::FrmOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void frm_loop(field::FrmOnv &bra, field::FrmOnv &ket, const PairBase::function_t &fn) override {
                base_t::loop_fn(bra, ket, fn);
            }
        };
    }

    namespace bos {
        struct Base : mbf_foreach::Base {
            const sys::bos::Sector m_sector;
            Base(const sys::bos::Sector& sector, size_t niter);     
        protected:
            void verify_mbf(const field::BosOnv& mbf) {
                REQUIRE_TRUE(mbf.m_basis==m_sector.m_basis, "given MBF has incorrect basis");
            }
        };

        struct GeneralClosed : Base {
            typedef basic_foreach::rtnd::Ordered<false, true> foreach_t;
            foreach_t m_foreach;

            GeneralClosed(const sys::bos::Sector& sector);

        public:
            template<typename fn_t>
            void loop_fn(field::BosOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                verify_mbf(mbf);
                auto basic_fn = [&mbf, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    mbf.set_ops(inds);
                    fn();
                };
                m_foreach.loop(basic_fn);
            }

        protected:
            void bos_loop(field::BosOnv &mbf, const function_t &fn) override;
        };

        struct GeneralOpen : Base {
            typedef basic_foreach::rtnd::Unrestricted foreach_t;
            foreach_t m_foreach;

            GeneralOpen(const sys::bos::Sector& sector);

        public:
            template<typename fn_t>
            void loop_fn(field::BosOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                verify_mbf(mbf);
                auto basic_fn = [&mbf, &fn](const basic_foreach::rtnd::inds_t &inds) {
                    mbf = inds;
                    fn();
                };
                m_foreach.loop(basic_fn);
            }

        protected:
            void bos_loop(field::BosOnv &mbf, const function_t &fn) override;

        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::BosOnv, foreach_t> {
            static_assert(std::is_base_of<bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::bos::Base");
            typedef mbf_foreach::Pair<field::BosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void bos_loop(field::BosOnv &bra, field::BosOnv &ket, const PairBase::function_t &fn) override {
                base_t::loop_fn(bra, ket, fn);
            }
        };
    }

    namespace frm_bos {
        struct Base : mbf_foreach::Base {
            Base(const sys::Sector& sector, size_t niter);
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
            struct Bases {
                const frm::Base & m_frm;
                const bos::Base & m_bos;
                Bases (const frm_foreach_t &frm, const bos_foreach_t &bos): m_frm(frm), m_bos(bos){}
            };

            static sys::Sector make_sector(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) {
                Bases bases(frm_foreach, bos_foreach);
                return {bases.m_frm.m_sector, bases.m_bos.m_sector};
            }

            static size_t make_niter(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) {
                Bases bases(frm_foreach, bos_foreach);
                return bases.m_frm.m_niter * bases.m_bos.m_niter;
            }

        public:
            Product(const frm_foreach_t &frm_foreach, const bos_foreach_t &bos_foreach) :
                    Base(make_sector(frm_foreach, bos_foreach), make_niter(frm_foreach, bos_foreach)),
                    m_frm_foreach(frm_foreach), m_bos_foreach(bos_foreach) {}


            template<typename fn_t>
            void loop_fn(field::FrmBosOnv &mbf, const fn_t &fn) {
                functor_utils::assert_prototype<prototype_t>(fn);
                auto frm_loop_fn = [this, &mbf, &fn]() {
                    auto bos_loop_fn = [&mbf, &fn]() {
                        fn();
                    };
                    m_bos_foreach.loop(mbf.m_bos, bos_loop_fn);
                };
                m_frm_foreach.loop(mbf.m_frm, frm_loop_fn);
            }

        protected:
            void frmbos_loop(field::FrmBosOnv &mbf, const function_t &fn) override {
                loop_fn(mbf, fn);
            }
        };

        /**
         * convenient partial specialization for fermions coupled to a closed bosonic quantum system
         * @tparam frm_foreach_t
         *  fermion foreach iterator type
         */
        template<typename frm_foreach_t>
        struct ClosedProduct : Product<frm_foreach_t, bos::GeneralClosed> {
            ClosedProduct(const frm_foreach_t &frm_foreach, const sys::bos::Sector& sector) :
                    Product<frm_foreach_t, bos::GeneralClosed>(frm_foreach, sector) {}
        };

        /**
         * convenient partial specialization for fermions coupled to an open bosonic quantum system
         * @tparam frm_foreach_t
         *  fermion foreach iterator type
         */
        template<typename frm_foreach_t>
        struct OpenProduct : Product<frm_foreach_t, bos::GeneralOpen> {
            OpenProduct(const frm_foreach_t &frm_foreach, const sys::bos::Sector& sector) :
                    Product<frm_foreach_t, bos::GeneralOpen>(frm_foreach, sector) {}
        };

        template<typename foreach_t>
        struct Pair : mbf_foreach::Pair<field::FrmBosOnv, foreach_t> {
            static_assert(std::is_base_of<frm_bos::Base, foreach_t>::value,
                          "template arg must be derived from mbf_foreach::frm_bos::Base");
            typedef mbf_foreach::Pair<field::FrmBosOnv, foreach_t> base_t;

            Pair(const foreach_t &foreach) : base_t(foreach) {}

        protected:
            void frmbos_loop(field::FrmBosOnv &bra, field::FrmBosOnv &ket, const PairBase::function_t &fn) override {
                base_t::loop_fn(bra, ket, fn);
            }
        };
    }
}


#endif //M7_MBFFOREACH_H
