//
// Created by rja on 28/03/2022.
//

#ifndef M7_CONNFOREACH_H
#define M7_CONNFOREACH_H

#include "BasicForeach.h"
#include "M7_lib/basis/Suites.h"

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


    namespace frm {
        struct Base : conn_foreach::Base {
            Base(size_t exsig, size_t nsite) :
                    conn_foreach::Base(exsig, {nsite, 0ul}) {
                REQUIRE_TRUE(exsig_utils::is_pure_frm(exsig), "excitation signature has boson operators");
            }
        };

        template<size_t nop>
        struct General : Base {
            General(size_t nsite) : Base(exsig_utils::encode(nop, nop, 0, 0), nsite) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                const auto& occs = src.m_decoded.m_simple_occs.get();
                const auto& vacs = src.m_decoded.m_simple_vacs.get();
                auto ann_fn = [&conn, &occs, &vacs, &fn](const ctnd::inds_t<nop>& ann_ops){
                    conn.m_ann.clear();
                    for (size_t iop=0ul; iop<nop; ++iop) conn.m_ann.add(occs[ann_ops[iop]]);
                    auto cre_fn = [&conn, &vacs, &fn](const ctnd::inds_t<nop>& cre_ops) {
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
        };

        template<size_t nop>
        struct Ms2Conserve : Base {
            Ms2Conserve(size_t nsite):
                Base(exsig_utils::encode(nop, nop, 0, 0), nsite) {}

        private:

            template<size_t nbeta, typename fn_t>
            void loop_one_beta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                const auto& occs = src.m_decoded.m_spin_occs.get();
                const auto& vacs = src.m_decoded.m_spin_vacs.get();
                constexpr auto nalpha = nop - nbeta;

                auto ann_alpha_fn = [&](const ctnd::inds_t<nalpha>& ann_alpha_ops) {
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

            template<size_t nbeta, typename fn_t>
            void loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn, tags::Bool<false>) {
                loop_one_beta_fn<nbeta, fn_t>(conn, src, fn);
                constexpr bool exit = nbeta == nop;
                loop_all_nbeta_fn<nbeta + 1, fn_t>(conn, src, fn, tags::Bool<exit>());
            }

            template<size_t nbeta, typename fn_t>
            void loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn, tags::Bool<true>) {}

        public:

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                /*
                 * to conserve 2*Ms, the number of beta electrons annihilated should equal the number of beta
                 * electrons created, so we need an outer loop over all numbers of betas, this is implemented in
                 * compile-time recursion
                 */
                loop_all_nbeta_fn<0, fn_t>(conn, src, fn, tags::Bool<false>());
            }

            template<typename fn_t>
            void loop_fn(const field::FrmOnv &src, const fn_t &fn) {
                loop_fn(m_conns.m_frmonv, src, fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t <conn::FrmOnv> &fn) override {
                Base::frm_loop(conn, src, fn);
            }
        };
    }

    namespace bos {
        struct Base : conn_foreach::Base {
            Base(size_t exsig, size_t nmode) :
                    conn_foreach::Base(exsig, {0ul, nmode}) {
                REQUIRE_TRUE(exsig_utils::is_pure_bos(exsig), "excitation signature has fermion operators");
            }
        };

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
};


#endif //M7_CONNFOREACH_H
