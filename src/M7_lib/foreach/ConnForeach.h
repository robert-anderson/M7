//
// Created by rja on 28/03/2022.
//

#ifndef M7_CONNFOREACH_H
#define M7_CONNFOREACH_H

#include "BasicForeach.h"
#include "M7_lib/basis/Suites.h"
#include "M7_lib/basis/Lattice.h"

namespace conn_foreach {
    using namespace basic_foreach;

    struct Base {
        const size_t m_exsig;
        const BasisData m_bd;
        suite::Conns m_conns;

        Base(size_t exsig, const BasisData &bd);

        virtual ~Base() {}

        template<typename conn_t>
        using function_t = std::function<void(const conn_t &)>;
    protected:
        virtual void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t<conn::FrmOnv> &fn) {};

        virtual void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t<conn::BosOnv> &fn) {};

        virtual void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src,
                                 const function_t<conn::FrmBosOnv> &fn) {};


    public:
        void loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t<conn::FrmOnv> &fn);

        void loop(const field::FrmOnv &src, const function_t<conn::FrmOnv> &fn);

        void loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t<conn::BosOnv> &fn);

        void loop(const field::BosOnv &src, const function_t<conn::BosOnv> &fn);

        void loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t<conn::FrmBosOnv> &fn);

        void loop(const field::FrmBosOnv &src, const function_t<conn::FrmBosOnv> &fn);
    };


    typedef std::unique_ptr<Base> base_ptr_t;
    typedef std::forward_list<base_ptr_t> base_list_t;

    namespace frm {
        struct Base : conn_foreach::Base {
            Base(size_t exsig, const FrmBasisData& bd) :
                    conn_foreach::Base(exsig, {bd, {}}) {
                REQUIRE_TRUE(exsig_utils::is_pure_frm(exsig), "excitation signature has boson operators");
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t <conn::FrmBosOnv> &fn) override;
        };

        template<size_t nop>
        struct General : Base {
            General(const FrmBasisData& bd) :
                Base(exsig_utils::encode(nop, nop, 0, 0), bd) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                const auto &occs = src.m_decoded.m_simple_occs.get();
                const auto &vacs = src.m_decoded.m_simple_vacs.get();
                auto ann_fn = [&conn, &occs, &vacs, &fn](const ctnd::inds_t<nop> &ann_ops) {
                    conn.m_ann.clear();
                    for (size_t iop = 0ul; iop < nop; ++iop) conn.m_ann.add(occs[ann_ops[iop]]);
                    auto cre_fn = [&conn, &vacs, &fn](const ctnd::inds_t<nop> &cre_ops) {
                        conn.m_cre.clear();
                        for (size_t iop = 0ul; iop < nop; ++iop) conn.m_cre.add(vacs[cre_ops[iop]]);
                        fn(conn);
                    };
                    ctnd::Ordered<nop, true, true> cre_foreach(vacs.size());
                    cre_foreach.loop(cre_fn);
                };
                ctnd::Ordered<nop, true, true> ann_foreach(occs.size());
                ann_foreach.loop(ann_fn);
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };

        template<size_t nop>
        struct Ms2Conserve : Base {
            Ms2Conserve(const FrmBasisData& bd):
                    Base(exsig_utils::encode(nop, nop, 0, 0), bd) {}

        private:

            template<typename fn_t, size_t nbeta, size_t nalpha>
            void loop_one_beta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn,
                                  tags::Int<nbeta>, tags::Int<nalpha>) {
                const auto &occs = src.m_decoded.m_spin_occs.get();
                const auto &vacs = src.m_decoded.m_spin_vacs.get();

                auto ann_alpha_fn = [&](const ctnd::inds_t<nalpha> &ann_alpha_ops) {
                    auto ann_beta_fn = [&](const ctnd::inds_t<nbeta> &ann_beta_ops) {
                        conn.m_ann.clear();
                        for (size_t iop = 0ul; iop < nalpha; ++iop) conn.m_ann.add(occs[0][ann_alpha_ops[iop]]);
                        for (size_t iop = 0ul; iop < nbeta; ++iop) conn.m_ann.add(occs[1][ann_beta_ops[iop]]);
                        auto cre_alpha_fn = [&](const ctnd::inds_t<nalpha> &cre_alpha_ops) {
                            auto cre_beta_fn = [&](const ctnd::inds_t<nbeta> &cre_beta_ops) {
                                conn.m_cre.clear();
                                for (size_t iop = 0ul; iop < nalpha; ++iop) conn.m_cre.add(vacs[0][cre_alpha_ops[iop]]);
                                for (size_t iop = 0ul; iop < nbeta; ++iop) conn.m_cre.add(vacs[1][cre_beta_ops[iop]]);
                                fn(conn);
                            };
                            ctnd::Ordered<nbeta, true, true> cre_beta_foreach(vacs.size(1));
                            nbeta ? cre_beta_foreach.loop(cre_beta_fn) : cre_beta_fn({});
                        };
                        ctnd::Ordered<nalpha, true, true> cre_alpha_foreach(vacs.size(0));
                        nalpha ? cre_alpha_foreach.loop(cre_alpha_fn) : cre_alpha_fn({});
                    };
                    ctnd::Ordered<nbeta, true, true> ann_beta_foreach(occs.size(1));
                    nbeta ? ann_beta_foreach.loop(ann_beta_fn) : ann_beta_fn({});
                };
                ctnd::Ordered<nalpha, true, true> ann_alpha_foreach(occs.size(0));
                // all loop bodies are explicitly called when the number of indices is zero
                nalpha ? ann_alpha_foreach.loop(ann_alpha_fn) : ann_alpha_fn({});
            }

            template<typename fn_t, size_t nbeta>
            void loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn, tags::Int<nbeta> tag) {
                static_assert(nop >= nbeta, "number of beta-spin operators cannot exceed excit level");
                loop_one_beta_fn<fn_t>(conn, src, fn, tag, tags::Int<nop - nbeta>());
                loop_all_nbeta_fn<fn_t>(conn, src, fn, tags::Int<nbeta + 1>());
            }

            template<typename fn_t>
            void
            loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn, tags::Int<nop + 1> tag) {}

        public:

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                /*
                 * to conserve 2*Ms, the number of beta electrons annihilated should equal the number of beta
                 * electrons created, so we need an outer loop over all numbers of betas, this is implemented in
                 * compile-time recursion
                 */
                loop_all_nbeta_fn<fn_t>(conn, src, fn, tags::Int<0>());
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };


        struct Hubbard : Base {
            const Lattice &m_lattice;

            Hubbard(const Lattice &lattice) : Base(exsig_utils::ex_single, {lattice.nsite()}),
                m_lattice(lattice) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                const auto &occs = src.m_decoded.m_simple_occs.get();
                for (const auto &occ: occs) {
                    conn.m_ann.clear();
                    conn.m_ann.add(occ);
                    auto ispin_occ = src.ispin(occ);
                    auto isite_occ = src.isite(occ);
                    auto coordinated_sites = m_lattice.m_sparse[isite_occ].first;
                    for (const auto &i: coordinated_sites) {
                        if (src.get({ispin_occ, i})) continue;
                        conn.m_cre.clear();
                        conn.m_cre.add({ispin_occ, i});
                        fn(conn);
                    }
                }
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };


        struct Heisenberg : Base {
            const Lattice &m_lattice;

            Heisenberg(const Lattice &lattice) :
                    Base(exsig_utils::ex_double, {lattice.nsite()}), m_lattice(lattice) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                // TODO: rewrite with m_decoded.m_alpha_only_occs when implemented
                const auto &occs = src.m_decoded.m_simple_occs.get();
                for (const auto &occ: occs) {
                    auto ispin_occ = src.ispin(occ);
                    if (ispin_occ) return; // all alpha bits have been dealt with
                    auto isite_occ = src.isite(occ);
                    // cannot exchange if the site is doubly occupied:
                    if (src.get({!ispin_occ, isite_occ})) continue;
                    auto coordinated_sites = m_lattice.m_sparse[isite_occ].first;
                    for (const auto &i: coordinated_sites) {
                        if (src.get({ispin_occ, i})) continue;
                        if (!src.get({!ispin_occ, i})) continue;
                        conn.m_ann.set({0, isite_occ}, {1, i});
                        conn.m_cre.set({0, i}, {1, isite_occ});
                        DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "incorrect excitation level");
                        fn(conn);
                    }
                }
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };
    }

    namespace bos {
        struct Base : conn_foreach::Base {
            Base(size_t exsig, const BosBasisData& bd) :
                    conn_foreach::Base(exsig, {{}, bd}){
                REQUIRE_TRUE(exsig_utils::is_pure_bos(exsig), "excitation signature has fermion operators");
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src,
                             const function_t <conn::FrmBosOnv> &fn) override {
                auto bos_fn = [&conn, &fn](const conn::BosOnv &bos_conn) {
                    fn(conn);
                };
                bos_loop(conn.m_bos, src.m_bos, bos_fn);

            }
        };

        struct Ann : Base {
            Ann(const BosBasisData& bd) : Base(exsig_utils::ex_0001, bd) {}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                conn.clear();
                const auto &occs = src.m_decoded.m_occ_modes.get();
                for (auto &imode: occs) {
                    conn.m_ann.set(imode);
                    fn(conn);
                }
            }

            template<typename fn_t>
            void loop_fn(const field::BosOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_bosonv, src, fn);
            }

        protected:
            void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t <conn::BosOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };

        struct Cre : Base {
            Cre(const BosBasisData& bd) : Base(exsig_utils::ex_0010, bd) {}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                conn.clear();
                for (size_t imode = 0ul; imode < src.m_bd.m_nmode; ++imode) {
                    if (size_t(src[imode] + 1) > m_bd.m_bos.m_nboson_max) continue;
                    conn.m_cre.set(imode);
                    fn(conn);
                }
            }

            template<typename fn_t>
            void loop_fn(const field::BosOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_bosonv, src, fn);
            }

        protected:
            void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t <conn::BosOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };

#if 0
        template<size_t nop_cre, size_t nop_ann>
        struct General : Base {

            General(size_t nmode, size_t nboson_max):
                Base(exsig_utils::encode(0, 0, nop_cre, nop_ann), nmode){}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                const auto& occs = src.m_decoded.m_occ_modes.get();

                auto nboson_max = m_nboson_max;
                auto ann_fn = [&](const ctnd::inds_t<nop_ann>& occ_inds) {
                    conn.m_ann.clear();
                    for (auto &ind: occ_inds) {
                        auto imode = occs[ind];
                        auto nop = conn.m_ann.add_to(imode);
                        // if the number of annihilations on the mode exceeds the occupation of the mode, then return
                        if (nop > src[imode]) return;
                    }
                    auto cre_fn = [&](const ctnd::inds_t<nop_ann> &mode_inds) {
                        conn.m_cre.clear();
                        for (auto &ind: mode_inds) {
                            auto imode = occs[ind];
                            auto nop = conn.m_ann.add_to(imode);
                            // if the number of annihilations on the mode exceeds the occupation of the mode, then return
                            if (nop + src[imode] > nboson_max) return;
                        }
                    };
                };
                ctnd::Ordered<nop_ann, false, true> ann_foreach(occs.size());
            }

            template<typename fn_t>
            void loop_fn(const field::BosOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_bosonv, src, fn);
            }

        protected:
            void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t<conn::BosOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };

#endif


#if 0
        template<size_t nop>
        struct GeneralClosed : Base {
            GeneralClosed(size_t nmode) : Base(exsig_utils::encode(0, 0, nop, nop), nmode) {}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                const auto& occs = src.m_decoded.m_occ_modes.get();
                auto ann_fn = [&conn, &occs, &fn](const ctnd::inds_t<nop>& ann_ops){
                    conn.m_ann.clear();
                    for (size_t iop=0ul; iop<nop; ++iop) conn.m_ann.add(occs[ann_ops[iop]]);
                    auto cre_fn = [&conn, &fn](const ctnd::inds_t<nop>& cre_ops) {
                        conn.m_cre.clear();
                        for (size_t iop=0ul; iop<nop; ++iop) conn.m_cre.add(vacs[cre_ops[iop]]);
                        fn(conn);
                    };
                    ctnd::Ordered<nop, true, true> cre_foreach(vacs.size());
                    cre_foreach.loop(cre_fn);
                };
                ctnd::Ordered<nop, true, true> ann_foreach(occs.size());
                ann_foreach.loop(ann_fn);
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                Base::frm_loop(conn, src, fn);
            }
        }
#endif
    }

    namespace frmbos {
        struct Base : conn_foreach::Base {
            const BasisData m_bd;

            Base(size_t exsig, const BasisData& bd) :
                    conn_foreach::Base(exsig, bd), m_bd(bd) {
                REQUIRE_TRUE(exsig_utils::decode_nfrm(exsig) && exsig_utils::decode_nbos(exsig),
                             "excitation signature is not that of a fermion-boson product");
            }
        };

        template<typename frm_t, typename bos_t>
        struct Product : Base {
            static_assert(std::is_base_of<frm::Base, frm_t>::value,
                          "template arg must be derived from conn_foreach::frm::Base");
            static_assert(std::is_base_of<bos::Base, bos_t>::value,
                          "template arg must be derived from conn_foreach::bos::Base");
        private:
            frm_t m_frm_foreach;
            bos_t m_bos_foreach;

            static const frm::Base& as_base(const frm_t& frm){
                return static_cast<const frm::Base&>(frm);
            }
            static const bos::Base& as_base(const bos_t& bos){
                return static_cast<const bos::Base&>(bos);
            }
            static size_t combined_exsig(const frm_t& frm, const bos_t& bos) {
                const frm::Base& frm_base = as_base(frm);
                auto nfrm_cre = exsig_utils::decode_nfrm_cre(frm_base.m_exsig);
                auto nfrm_ann = exsig_utils::decode_nfrm_ann(frm_base.m_exsig);
                const bos::Base& bos_base = as_base(bos);
                auto nbos_cre = exsig_utils::decode_nbos_cre(bos_base.m_exsig);
                auto nbos_ann = exsig_utils::decode_nbos_ann(bos_base.m_exsig);
                return exsig_utils::encode(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann);
            }
            static BasisData combined_bd(const frm_t& frm, const bos_t& bos) {
                const frm::Base& frm_base = as_base(frm);
                const bos::Base& bos_base = as_base(bos);
                return {frm_base.m_bd.m_frm, bos_base.m_bd.m_bos};
            }

        public:
            Product(frm_t &&frm_foreach, bos_t &&bos_foreach) :
                    Base(combined_exsig(frm_foreach, bos_foreach), combined_bd(frm_foreach, bos_foreach)),
                    m_frm_foreach(std::move(frm_foreach)), m_bos_foreach(std::move(bos_foreach)) {
            }

            template<typename fn_t>
            void loop_fn(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const fn_t &fn) {
                conn.clear();
                auto frm_fn = [&](const conn::FrmOnv& frm_conn){
                    auto bos_fn = [&](const conn::BosOnv& bos_conn) {
                        fn(conn);
                    };
                    m_bos_foreach.loop_fn(conn.m_bos, src.m_bos, bos_fn);
                };
                m_frm_foreach.loop_fn(conn.m_frm, src.m_frm, frm_fn);
            }

            template<typename fn_t>
            void loop_fn(const field::FrmBosOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmbosonv, src, fn);
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src,
                             const function_t<conn::FrmBosOnv> &fn) override {
                loop_fn(conn, src, fn);
            }
        };
    }
}

#endif //M7_CONNFOREACH_H
