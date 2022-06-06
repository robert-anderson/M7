//
// Created by Robert J. Anderson on 12/8/21.
//

#ifndef M7_HUBBARDFRMHAM_H
#define M7_HUBBARDFRMHAM_H

#include "M7_lib/foreach/Foreach.h"
#include "M7_lib/nd/NdFormatD.h"
#include "M7_lib/basis/Lattice.h"
#include "M7_lib/linalg/Dense.h"
#include "M7_lib/foreach/ConnForeach.h"

#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct HubbardFrmHam : FrmHam {
    /**
     * on-site repulsion scalar in units of the hopping
     */
    const defs::ham_t m_u;

    HubbardFrmHam(defs::ham_t u, const std::shared_ptr<lattice::Base>& lattice);

    HubbardFrmHam(opt_pair_t opts);

    defs::ham_t get_coeff_1100(size_t a, size_t i) const override;

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    void log_data() const override;

    excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& opts) const override;

    conn_foreach_list_t make_foreach_iters() const override;
};


#endif //M7_HUBBARDFRMHAM_H
