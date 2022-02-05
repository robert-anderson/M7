//
// Created by anderson on 12/8/21.
//

#ifndef M7_HUBBARDFRMHAM_H
#define M7_HUBBARDFRMHAM_H

#include <src/core/util/Foreach.h>
#include <src/core/nd/NdFormatD.h>
#include <src/core/basis/Lattice.h>
#include "FrmHam.h"
#include "src/core/linalg/Dense.h"

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

    defs::ham_t get_coeff_1100(const size_t &i, const size_t &j) const override;

    defs::ham_t get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    void log_data() const override;
};


#endif //M7_HUBBARDFRMHAM_H
