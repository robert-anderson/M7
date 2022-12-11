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
        uint_t get_occ_uniform(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const;

        prob_t prob_uniform(const field::FrmOnv &src, const conn::FrmOnv &conn) const;

        virtual uint_t get_occ(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const = 0;

    public:

        bool draw_frm(OpSig exsig, const field::FrmOnv &src, prob_t &prob, conn::FrmOnv &conn) override;

        uint_t approx_nconn(OpSig exsig, sys::Particles particles) const override;

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

        HubbardPreferDoubleOcc(const FrmHam &h, PRNG &prng, prob_t doub_occ_u_fac);
    protected:

        prob_t combined_occ_prob(uint_t nelec, uint_t nelec_doub_occ) const;

        uint_t get_occ(uint32_t rand, const field::FrmOnv &src, prob_t &prob) const override;

    public:
        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;
    };
}

#endif //M7_EXGEN_HUBBARD_H
