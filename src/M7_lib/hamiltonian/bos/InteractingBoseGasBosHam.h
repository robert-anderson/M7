//
// Created by Robert J. Anderson on 07/03/2022.
//

#ifndef M7_INTERACTINGBOSEGASBOSHAM_H
#define M7_INTERACTINGBOSEGASBOSHAM_H

#include "M7_lib/conf/HamiltonianConf.h"
#include "M7_lib/basis/Planewaves.h"
#include "M7_lib/nd/NdFormatD.h"

#include "M7_lib/hamiltonian/bos/BosHam.h"

struct InteractingBoseGasBosHam : BosHam {
    const Planewaves m_planewaves;
    const ham_comp_t m_ek_scale;

    InteractingBoseGasBosHam(uint_t ndim, uint_t nwave, ham_comp_t ek_scale);

    InteractingBoseGasBosHam(init_opts_t opts);

    ham_t get_coeff_0011(uint_t /*i*/, uint_t /*j*/) const override {
        return 0.0;
    }

    ham_t get_coeff_0022(uint_t i, uint_t j, uint_t k, uint_t l) const override;

    ham_t get_element_0000(const field::BosOnv& onv) const override;

    ham_t get_element_0011(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const override {
        return 0.0;
    }

    ham_t get_element_0022(const field::BosOnv& /*onv*/, const conn::BosOnv& conn) const override;
};


#endif //M7_INTERACTINGBOSEGASBOSHAM_H
