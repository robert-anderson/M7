//
// Created by anderson on 12/8/21.
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
    /**
     * encodes all connection information between the flattened indices of the Hubbard lattice sites
     */
    Lattice m_lattice;
    const NdFormatD& m_format;
    const std::vector<int>& m_bcs;
    /**
     * whether the model meets the sign problem-free conditions
     */
    const bool m_spf;

private:

    /**
     * determine whether this hubbard hamiltonian is sign-problematic
     * @return
     *  true if the model is sign problem-free
     */
    bool sign_problem() const;


public:

    HubbardFrmHam(defs::ham_t u, Lattice lattice, int ms2_restrict, int charge);

    HubbardFrmHam(const fciqmc_config::FermionHamiltonian &opts);

    defs::ham_t get_coeff_1100(size_t a, size_t i) const override;

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    void log_data() const override;

    excit_gen_list_t make_excit_gens(PRNG& prng, const fciqmc_config::Propagator& opts) override;

    conn_iter_ptr_list_t make_conn_iters() override;
};


#endif //M7_HUBBARDFRMHAM_H
