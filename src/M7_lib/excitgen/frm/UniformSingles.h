//
// Created by Robert J. Anderson on 04/04/2022.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FrmExcitGen.h"

struct UniformSingles : FrmExcitGen {
    UniformSingles(const FrmHam &h, PRNG &prng):
            FrmExcitGen(h, prng, {exsig_utils::ex_single}, "uniform"){}

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) override;

    static bool draw_spin_conserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

    static bool draw_spin_nonconserve_fn(PRNG &prng, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn);

private:
    defs::prob_t prob_spin_conserve_fn(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
        const auto &nonempty_pairs = src.m_decoded.m_nonempty_pair_labels.get();
        if (nonempty_pairs.empty()) return 0.0;
        const auto label = src.m_decoded.m_spin_sym_occs.label(conn.m_cre[0]);
        const auto &occs = src.m_decoded.m_spin_sym_occs.get()[label];
        const auto &vacs = src.m_decoded.m_spin_sym_vacs.get()[label];
        return 1.0 / (nonempty_pairs.size() * occs.size() * vacs.size());
    }

    defs::prob_t prob_spin_nonconserve_fn(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
        const auto nocc = src.m_decoded.m_simple_occs.get().size();
        const auto nvac = src.m_basis.m_nspinorb - nocc;
        return 1.0/(nocc*nvac);
    }

public:

    defs::prob_t prob_frm(const field::FrmOnv &src, const conn::FrmOnv &conn) const override;

    size_t approx_nconn(size_t exsig, sys::Particles particles) const override;
};


#endif //M7_UNIFORMSINGLES_H
