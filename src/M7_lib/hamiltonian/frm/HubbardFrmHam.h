//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_HUBBARDFRMHAM_H
#define M7_HUBBARDFRMHAM_H

#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct HubbardFrmHam : FrmHam {
    /**
     * on-site repulsion scalar in units of the hopping
     */
    const ham_comp_t m_u;

    HubbardFrmHam(ham_comp_t u, const std::shared_ptr<lattice::Lattice>& lattice);

    HubbardFrmHam(opt_pair_t opts);

    ham_t get_coeff_1100(uint_t a, uint_t i) const override;

    ham_t get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const override;

    ham_t get_element_0000(const field::FrmOnv& onv) const override;

    ham_t get_element_1100(const field::FrmOnv& onv, const conn::FrmOnv& conn) const override;

    ham_t get_element_2200(const field::FrmOnv& /*onv*/, const conn::FrmOnv& /*conn*/) const override {
        return 0;
    }

    void log_data() const override;

    excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& opts) const override;

    conn_foreach_list_t make_foreach_iters() const override;
};


#endif //M7_HUBBARDFRMHAM_H
