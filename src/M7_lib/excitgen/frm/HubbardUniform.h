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
        virtual void set_non_uniform_occ(uint_t &occ, const field::FrmOnv &src, prob_t &prob) const = 0;

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
        void set_non_uniform_occ(uint_t& /*occ*/, const field::FrmOnv& src, prob_t &prob) const override;

    public:
        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;
    };
    /**
     * biases the creation operator selection so as to prefer electrons which are on doubly-occupied sites
     */
//    struct HubbardPreferDoubleOcc : HubbardBase {
//
//        HubbardPreferDoubleOcc(const FrmHam &h, PRNG &prng):
//            HubbardBase(h, prng, "doubly occupied site-preferring hubbard hopping"){}
//    protected:
//        void set_non_uniform_occ(uint_t& /*occ*/, const field::FrmOnv& src, prob_t &prob) const override {
//            prob = 1/double(src.m_decoded.m_simple_occs.get().size());
//        }
//
//    public:
//        prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;
//    };
}

#endif //M7_EXGEN_HUBBARD_H
