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
    const defs::ham_t m_ek_scale;

    InteractingBoseGasBosHam(uint_t ndim, uint_t nwave, defs::ham_t ek_scale):
        BosHam(Planewaves::size(ndim, nwave)), m_planewaves(ndim, nwave), m_ek_scale(ek_scale){}

    InteractingBoseGasBosHam(opt_pair_t opts):
            InteractingBoseGasBosHam(
                    opts.m_ham.m_interacting_bose_gas.m_ndim,
                    opts.m_ham.m_interacting_bose_gas.m_nwave,
                    opts.m_ham.m_interacting_bose_gas.m_ek_scale){}

    defs::ham_t get_coeff_0011(uint_t /*i*/, uint_t /*j*/) const override {
        return 0.0;
    }

    defs::ham_t get_coeff_0022(uint_t i, uint_t j, uint_t k, uint_t l) const override {
        return m_planewaves.conserving(i, j, k, l) ? 1.0 : 0.0;
    }

    defs::ham_t get_element_0000(const field::BosOnv& onv) const override {
        // total linear momentum
        defs::ham_t tot = 0.0;
        for (uint_t imode=0ul; imode<m_basis.m_nmode; ++imode){
            if (!onv[imode]) continue;
            tot+=m_planewaves.kinetic_energy(imode)*onv[imode];
            //tot+=0.5*onv[imode]*(onv[imode]-1);
        }
        return tot * m_ek_scale;
    }

    defs::ham_t get_element_0011(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const override {
        return 0.0;
    }

    defs::ham_t get_element_0022(const field::BosOnv& /*onv*/, const conn::BosOnv& conn) const override {
        const auto i = conn.m_cre.get_imode(0);
        const auto j = conn.m_cre.get_imode(1);
        const auto k = conn.m_ann.get_imode(0);
        const auto l = conn.m_ann.get_imode(1);
        return get_coeff_0022(i, j, k, l);
    }
};


#endif //M7_INTERACTINGBOSEGASBOSHAM_H
