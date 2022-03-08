//
// Created by rja on 07/03/2022.
//

#ifndef M7_INTERACTINGBOSEGASBOSHAM_H
#define M7_INTERACTINGBOSEGASBOSHAM_H

#include <src/core/config/Hamiltonian.h>
#include <src/core/nd/NdFormatD.h>
#include "BosHam.h"

struct InteractingBoseGasBosHam : BosHam {
    const size_t m_nplanewave;
    /**
     * number of k-points per dimension
     */
    const size_t m_nkpoint;
    const NdFormatD m_format;
    const defs::ham_t m_ek_scale;

protected:
    /**
     * working arrays into which the mode indices are decoded
     */
    mutable std::vector<int> m_ikpoints_i_work, m_ikpoints_j_work, m_ikpoints_k_work, m_ikpoints_l_work;

public:
    void get_ikpoints(const size_t& imode, std::vector<int>& ikpoints) const {
        /*
         *  0  1  2  3  4  5  6
         * -3 -2 -1  0  1  2  3
         */
        m_format.decode_flat(imode, ikpoints);
        for (auto& ikpoint: ikpoints) ikpoint-=m_nplanewave;
    }

    defs::ham_t get_ksquare(const size_t& imode) const {
        m_format.decode_flat(imode, m_ikpoints_i_work);
        defs::ham_t tot = 0.0;
        for (const auto& ikpoint: m_ikpoints_i_work) {
            std::cout << ikpoint << std::endl;
            tot+=ikpoint*ikpoint;
        }
        return tot;
    }

    InteractingBoseGasBosHam(size_t nboson, size_t ndim, size_t nplanewave, defs::ham_t ek_scale):
            BosHam(std::pow(2*nplanewave+1ul, ndim), nboson),
            m_nplanewave(nplanewave), m_nkpoint(2*nplanewave+1), m_format(ndim, m_nkpoint), m_ek_scale(ek_scale){}

    InteractingBoseGasBosHam(const fciqmc_config::BosonHamiltonian &opts):
            InteractingBoseGasBosHam(opts.m_nboson, opts.m_interacting_bose_gas.m_ndim,
                 opts.m_interacting_bose_gas.m_nplanewave, opts.m_interacting_bose_gas.m_ek_scale){}


    defs::ham_t get_coeff_0011(const size_t &i, const size_t &j) const override {
        return 0.0;
    }

    defs::ham_t get_coeff_0022(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override {
        get_ikpoints(i, m_ikpoints_i_work);
        get_ikpoints(j, m_ikpoints_j_work);
        get_ikpoints(k, m_ikpoints_k_work);
        get_ikpoints(l, m_ikpoints_l_work);
        for (size_t idim=0ul; idim<m_format.m_nind; ++idim){
            auto cre_mom = m_ikpoints_i_work[idim]+m_ikpoints_j_work[idim];
            auto ann_mom = m_ikpoints_k_work[idim]+m_ikpoints_l_work[idim];
            // linear momentum must be conserved in all directions
            if (cre_mom!=ann_mom) return 0.0;
        }
        return 1.0;
    }

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override {
        // total linear momentum
        defs::ham_t tot = 0.0;
        for (size_t imode=0ul; imode<m_nmode; ++imode){
            if (!onv[imode]) continue;
            tot+=get_ksquare(imode)*onv[imode];
            std::cout << "" << std::endl;
            //tot+=0.5*onv[imode]*(onv[imode]-1);
        }
        return tot * m_ek_scale;
    }

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0.0;
    }

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        const auto i = conn.m_cre.get_imode(0);
        const auto j = conn.m_cre.get_imode(1);
        const auto k = conn.m_ann.get_imode(0);
        const auto l = conn.m_ann.get_imode(1);
        return get_coeff_0022(i, j, k, l);
    }
};


#endif //M7_INTERACTINGBOSEGASBOSHAM_H
