//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include <src/core/config/FciqmcConfig.h>
#include "Propagator.h"

class ExactPropagator : public Propagator {
    /**
     * sends the generated excitations even if the corresponding hamiltonian matrix element is zero. Useful for testing
     * rank-2 RDMs since these spawns will make the exact contributions.
     */
    const bool m_only_nonzero_h_spawns;

    defs::ham_t off_diagonal_bosons(const Hamiltonian<0> &ham, conn::FrmOnv &conn,
                                    const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv, const size_t &occ, int change){
        return 0;
    }

#if 0
    defs::ham_t off_diagonal_bosons(const Hamiltonian<1> &ham, conn::FrmBosOnv &conn,
                                    const fields::FrmBosOnv &src_onv, fields::FrmBosOnv &dst_onv, const size_t &occ, int change){
        const size_t imode = occ < ham.nsite() ? occ : occ-ham.nsite();
        if (src_onv.m_bos[imode] == 0 && (change < 0)) return 0.0;
        else if (src_onv.m_bos[imode] == ham.nboson_cutoff() && (change > 0)) return 0.0;

        conn.zero();
        conn.m_bonvconn.add(imode, change);
        dst_onv.zero();
        conn.apply(src_onv, dst_onv);
        ASSERT(src_onv.m_frm == dst_onv.m_frm);
        auto com = dst_onv.m_bos[imode];
        if (change<0) com+=change;
        auto helem = ham.bc().get_element_1(imode, imode, com);

        ASSERT(consts::floats_equal(helem, ham.bc().v(imode, imode, imode)*std::sqrt(com+1)))
        return helem;
    }
#endif


public:
    ExactPropagator(const Hamiltonian<> &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt,
                    bool only_nonzero_h_spawns=true);

    void diagonal(Wavefunction &wf, const size_t& ipart) override;

    void off_diagonal(Wavefunction &wf, const size_t& ipart) override;

    bool is_exact() const override;

};

#endif //M7_EXACTPROPAGATOR_H
