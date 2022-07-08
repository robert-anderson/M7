//
// Created by Robert J. Anderson on 04/04/2022.
//

#ifndef M7_EXGEN_HUBBARD_H
#define M7_EXGEN_HUBBARD_H

#include "FrmExcitGen.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

namespace exgen {
    /**
     * excitation generator for the N-dimensional Hubbard model
     */
    struct HubbardBase : FrmLatticeExcitGen {

        HubbardBase(const FrmHam &h, PRNG &prng, str_t description);

        virtual ~HubbardBase() {}

    protected:
        uint_t get_occ_uniform(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const {
            const auto &occs = src.m_decoded.m_simple_occs.get();
            const auto nconn_lcm = m_h.m_basis.m_lattice->m_lcm_le_nadj_max;
            const auto nelec = occs.size();
            prob = 1/double(nelec);
            return occs[rand / nconn_lcm];
        }

        prob_t prob_uniform(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
            const auto &occs = src.m_decoded.m_simple_occs.get();
            const auto nelec = occs.size();
            const auto isite = src.m_basis.isite(conn.m_ann[0]);
            const auto ispin = src.m_basis.ispin(conn.m_ann[0]);
            /*
             * fill the working vector with the adjacency of the occupied site
             */
            m_h.m_basis.m_lattice->get_adj_row(isite, m_work_adj_row);
            /*
             * when selecting the vacant site, skip the sites which are occupied in the same spin channel as the chosen electron
             */
            set_valid_adj_vacant(src, ispin);
            return 1.0/(m_valid_in_adj_row.size()*nelec);
        }

        virtual uint_t get_occ(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const = 0;

    public:

        bool draw_frm(uint_t exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) override;

        uint_t approx_nconn(uint_t exsig, sys::Particles particles) const override;

    };
    /**
     * only constrains the creation operator drawn based on the occupation of the adjacent sites of the chosen electron
     */
    struct HubbardUniform : HubbardBase {

        HubbardUniform(const FrmHam &h, PRNG &prng): HubbardBase(h, prng, "uniform hubbard hopping"){}
    protected:
        uint_t get_occ(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const override{
            return get_occ_uniform(rand, src, prob);
        }

    public:
        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override {
            return prob_uniform(src, conn);
        }
    };
    /**
     * biases the creation operator selection so as to prefer electrons which are on doubly-occupied sites
     */
    struct HubbardPreferDoubleOcc : HubbardBase {
        /**
         * in the event that there are any doubly occupied sites, decide whether to draw an electron from the set of
         * double occupations with this probability:
         * prob_doub_occ : prob_uniform = doub_occ_u_fac * U : 1
         * i.e. prob_doub_occ = 1/(1/(doub_occ_u_fac * U)+1)
         */
        const prob_t m_prob_doub_occ;

        HubbardPreferDoubleOcc(const FrmHam &h, PRNG &prng, prob_t doub_occ_u_fac):
            HubbardBase(h, prng, "doubly occupied site-preferring hubbard hopping"),
            m_prob_doub_occ(1.0/(1.0+1.0/(doub_occ_u_fac*m_h.as<HubbardFrmHam>()->m_u))){}
    protected:

        prob_t combined_occ_prob(uint_t nelec, uint_t nelec_doub_occ) const {
            // prob = prob_uni + prob_pref
            return (1.0 - m_prob_doub_occ)/nelec + m_prob_doub_occ/nelec_doub_occ;
        }

        uint_t get_occ(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const override {
            /*
             * total number of electrons
             */
            const auto nelec = src.m_decoded.m_simple_occs.get().size();
            const auto& doub_occs = src.m_decoded.m_doubly_occ_sites.get();
            /*
             * if there are no double occupations, then do uniform selection without alteration of prob
             */
            if (doub_occs.empty()) return get_occ_uniform(rand, src, prob);
            /*
             * number of electrons in doubly occupied sites
             */
            const auto nelec_doub_occ = 2*doub_occs.size();
            /*
             * it should still be possible to draw any electron
             */
            if (m_prng.draw_float() > m_prob_doub_occ) {
                /*
                 * uniform branch: draw using the uniform strategy
                 */
                const auto occ = get_occ_uniform(rand, src, prob);
                const auto isite = src.m_basis.isite(occ);
                const auto ispin = src.m_basis.ispin(occ);
                if (src.get({!ispin, isite})) {
                    DEBUG_ASSERT_EQ(src.site_nocc(isite), 2ul, "should have a double occupation");
                    /*
                     * this electron on a doubly-occupied site could have been drawn either by this way (uniform branch)
                     * or by the preferential branch
                     */
                    prob = combined_occ_prob(nelec, nelec_doub_occ);
                }
                else {
                    /*
                     * there is only one way to draw this electron in a singly-occupied site, but the probability must
                     * be scaled by the probability that the uniform draw was attempted in the presence of a double
                     * occupation
                     */
                    prob /= 1.0 - m_prob_doub_occ;
                }
                return occ;
            }
            /*
             * preferential branch
             */
            rand = m_prng.draw_uint(nelec_doub_occ);
            const size_t ispin = rand%2;
            const auto isite = doub_occs[rand/2];
            prob = combined_occ_prob(nelec, nelec_doub_occ);
            return src.m_basis.ispinorb(ispin, isite);
        }

    public:
        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override {
            const auto& doub_occs = src.m_decoded.m_doubly_occ_sites.get();
            if (doub_occs.empty()) return prob_uniform(src, conn);
            return combined_occ_prob(src.m_decoded.m_simple_occs.get().size(), doub_occs.size());
        }
    };
}

#endif //M7_EXGEN_HUBBARD_H
