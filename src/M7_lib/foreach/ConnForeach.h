//
// Created by Robert J. Anderson on 28/03/2022.
//

#ifndef M7_CONNFOREACH_H
#define M7_CONNFOREACH_H

#include "BasicForeach.h"
#include "M7_lib/basis/Suites.h"
#include "M7_lib/basis/Lattice.h"
#include "M7_lib/util/Exsig.h"

/**
 * A collection of iterators over connections with different constraints. These are divided into those which loop over
 * purely fermionic (frm), purely bosonic (bos), and fermion-boson product (frmbos).
 *
 * The foreach iterators defined here are run time polymorphic, and all reference to the MBF type is missing from the
 * base class so that the MBF type only needs to be specified in the calling scope when calling the overloaded loop
 * methods of the polymorphic conn_foreach::Base class.
 *
 * The derived classes are also compile time polymorphic through the definition of "loop_fn" templated methods, to which
 * the virtual methods should delegate. This means that when the calling scope has access to a specific derived type of
 * conn_foreach::Base, it can call the functor-parametrised loop_fn directly, which should exhibit performance benefits
 * over the runtime polymorphic use of "loop" (which delegates to the the virtual "frm_loop", "bos_loop", "frmbos_loop"
 * methods).
 */
namespace conn_foreach {
    using namespace basic_foreach;

    struct Base {
        const size_t m_exsig;

        Base(size_t exsig);

        virtual ~Base() {}


        using function_t = std::function<void()>;
    protected:
        virtual void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn) {};

        virtual void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t &fn) {};

        virtual void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t &fn) {};


    public:
        void loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn);

        void loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t &fn);

        void loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t &fn);
    };


    typedef std::unique_ptr<Base> base_ptr_t;
    typedef std::forward_list<base_ptr_t> base_list_t;

    namespace frm {
        struct Base : conn_foreach::Base {
            Base(size_t exsig) : conn_foreach::Base(exsig) {
                REQUIRE_TRUE(utils::exsig::is_pure_frm(exsig), "excitation signature has boson operators");
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t &fn) override;
        };

        template<size_t nop>
        struct General : Base {
            General() : Base(utils::exsig::encode(nop, nop, 0, 0)) {}

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
                        fn();
                    };
                    ctnd::Ordered<nop, true, true> cre_foreach(vacs.size());
                    cre_foreach.loop(cre_fn);
                };
                ctnd::Ordered<nop, true, true> ann_foreach(occs.size());
                ann_foreach.loop(ann_fn);
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };

        template<size_t nop>
        struct Ms2Conserve : Base {
            Ms2Conserve(): Base(utils::exsig::encode(nop, nop, 0, 0)) {}

        private:

            template<typename fn_t, size_t nbeta, size_t nalpha>
            void loop_one_beta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn,
                                  utils::tag::Int<nbeta>, utils::tag::Int<nalpha>) {
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
                                fn();
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
            void loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn, utils::tag::Int<nbeta> tag) {
                static_assert(nop >= nbeta, "number of beta-spin operators cannot exceed excit level");
                loop_one_beta_fn<fn_t>(conn, src, fn, tag, utils::tag::Int<nop - nbeta>());
                loop_all_nbeta_fn<fn_t>(conn, src, fn, utils::tag::Int<nbeta + 1>());
            }

            template<typename fn_t>
            void loop_all_nbeta_fn(conn::FrmOnv &conn, const field::FrmOnv &src,
                                   const fn_t &fn, utils::tag::Int<nop + 1> tag) {}

        public:

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                utils::functor::assert_prototype<void()>(fn);

                /*
                 * to conserve 2*Ms, the number of beta electrons annihilated should equal the number of beta
                 * electrons created, so we need an outer loop over all numbers of betas, this is implemented in
                 * compile-time recursion
                 */
                loop_all_nbeta_fn<fn_t>(conn, src, fn, utils::tag::Int<0>());
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };


        struct Hubbard : Base {
            Hubbard() : Base(utils::exsig::ex_single) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                utils::functor::assert_prototype<void()>(fn);
                const auto lattice = src.m_basis.m_lattice;
                REQUIRE_TRUE(lattice.get(), "Hubbard model requires the basis to have a lattice defined");
                const auto &occs = src.m_decoded.m_simple_occs.get();
                lattice::adj_row_t adj_row;
                for (const auto &occ: occs) {
                    conn.m_ann.clear();
                    conn.m_ann.add(occ);
                    auto ispin_occ = src.m_basis.ispin(occ);
                    auto isite_occ = src.m_basis.isite(occ);
                    lattice->get_adj_row(isite_occ, adj_row);
                    for (const auto &adj_elem: adj_row) {
                        auto vac = src.m_basis.ispinorb(ispin_occ, adj_elem.m_isite);
                        if (src.get(vac)) continue;
                        conn.m_cre.clear();
                        conn.m_cre.add(vac);
                        fn();
                    }
                }
            }

        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };


        struct Heisenberg : Base {
            Heisenberg() : Base(utils::exsig::ex_double) {}

            template<typename fn_t>
            void loop_fn(conn::FrmOnv &conn, const field::FrmOnv &src, const fn_t &fn) {
                utils::functor::assert_prototype<void()>(fn);
                // TODO: rewrite with m_decoded.m_alpha_only_occs when implemented
                const auto lattice = src.m_basis.m_lattice;
                REQUIRE_TRUE(lattice.get(), "Hubbard model requires the basis to have a lattice defined");
                const auto &occs = src.m_decoded.m_simple_occs.get();
                lattice::adj_row_t adj_row;
                for (const auto &occ: occs) {
                    auto ispin_occ = src.m_basis.ispin(occ);
                    if (ispin_occ) return; // all alpha bits have been dealt with
                    auto isite_occ = src.m_basis.isite(occ);
                    // cannot exchange if the site is doubly occupied:
                    if (src.get({!ispin_occ, isite_occ})) continue;
                    lattice->get_adj_row(isite_occ, adj_row);
                    for (const auto &adj_elem: adj_row) {
                        const auto i = adj_elem.m_isite;
                        if (src.get({ispin_occ, i})) continue;
                        if (!src.get({!ispin_occ, i})) continue;
                        conn.m_ann.set({0, isite_occ}, {1, i});
                        conn.m_cre.set({0, i}, {1, isite_occ});
                        DEBUG_ASSERT_EQ(conn.exsig(), utils::exsig::ex_double, "incorrect excitation level");
                        fn();
                    }
                }
            }


        protected:
            void frm_loop(conn::FrmOnv &conn, const field::FrmOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };
    }

    namespace bos {
        struct Base : conn_foreach::Base {
            Base(size_t exsig): conn_foreach::Base(exsig){
                REQUIRE_TRUE(utils::exsig::is_pure_bos(exsig), "excitation signature has fermion operators");
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t &fn) override {
                bos_loop(conn.m_bos, src.m_bos, fn);

            }
        };

        struct Ann : Base {
            Ann() : Base(utils::exsig::ex_0001) {}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                utils::functor::assert_prototype<void()>(fn);
                conn.clear();
                const auto &occs = src.m_decoded.m_occ_modes.get();
                for (auto &imode: occs) {
                    conn.m_ann.set(imode);
                    fn();
                }
            }

        protected:
            void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };

        struct Cre : Base {
            Cre(): Base(utils::exsig::ex_0010) {}

            template<typename fn_t>
            void loop_fn(conn::BosOnv &conn, const field::BosOnv &src, const fn_t &fn) {
                utils::functor::assert_prototype<void()>(fn);
                conn.clear();
                for (size_t imode = 0ul; imode < src.m_size; ++imode) {
                    if (size_t(src[imode] + 1) > src.m_basis.m_occ_cutoff) continue;
                    conn.m_cre.set(imode);
                    fn();
                }
            }

        protected:
            void bos_loop(conn::BosOnv &conn, const field::BosOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };
    }

    namespace frmbos {
        struct Base : conn_foreach::Base {
            Base(size_t exsig) : conn_foreach::Base(exsig) {
                REQUIRE_TRUE(utils::exsig::decode_nfrm(exsig) && utils::exsig::decode_nbos(exsig),
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

            static size_t combined_exsig() {
                const frm::Base frm = frm_t();
                auto nfrm_cre = utils::exsig::decode_nfrm_cre(frm.m_exsig);
                auto nfrm_ann = utils::exsig::decode_nfrm_ann(frm.m_exsig);
                const bos::Base bos = bos_t();
                auto nbos_cre = utils::exsig::decode_nbos_cre(bos.m_exsig);
                auto nbos_ann = utils::exsig::decode_nbos_ann(bos.m_exsig);
                return utils::exsig::encode(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann);
            }

        public:
            Product(): Base(combined_exsig()){}

            template<typename fn_t>
            void loop_fn(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const fn_t &fn) {
                conn.clear();
                auto frm_fn = [&](){
                    auto bos_fn = [&]() {
                        fn();
                    };
                    m_bos_foreach.loop_fn(conn.m_bos, src.m_bos, bos_fn);
                };
                m_frm_foreach.loop_fn(conn.m_frm, src.m_frm, frm_fn);
            }

        protected:
            void frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src, const function_t &fn) override {
                loop_fn(conn, src, fn);
            }
        };
    }
}

#endif //M7_CONNFOREACH_H
