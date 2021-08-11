//
// Created by rja on 11/08/2021.
//

#ifndef M7_BILINEARESTIMATORGROUP_H
#define M7_BILINEARESTIMATORGROUP_H

#include <utility>

#include "src/core/field/Fields.h"
#include "Rdm.h"
#include "SpectralMoment.h"


struct BilinearEstimatorGroup {

    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_rdm_exsigs;
    std::array<std::unique_ptr<SpectralMoment>, defs::nexsig> m_specmoms;
    const defs::inds m_specmom_exsigs;

    bool any_rdms() const {
        return !m_rdm_exsigs.empty();
    }

    bool any_specmoms() const {
        return !m_specmom_exsigs.empty();
    }

    /**
     * replication is needed if we are estimating any bilinears stochastically
     * @param stoch
     *  true if the wavefunction is propagated by a stochastic process
     * @return
     *  true if replication of the walker populations is a statistical necessity
     */
    bool need_replication(bool stoch) const {
        if (!stoch) return false;
        return any_rdms() || any_specmoms();
    }

    /**
     * parents only need to be sent in spawning table if we're going to be accumulating RDMs at some point
     */
    bool need_send_parents() const {
        return any_rdms();
    }

    BilinearEstimatorGroup(defs::inds rdm_exsigs, defs::inds specmom_exsigs) :
        m_rdm_exsigs(std::move(rdm_exsigs)), m_specmom_exsigs(std::move(specmom_exsigs)) {
        for (auto exsig: m_rdm_exsigs) {
            REQUIRE_TRUE(exsig, "multidimensional estimators require a nonzero number of SQ operator indices");
            if (conn_utils::pure_frm(exsig))
                m_rdms[exsig] = std::unique_ptr<Rdm>(new FrmRdm(exsig));
            else if (conn_utils::pure_bos(exsig))
                m_rdms[exsig] = std::unique_ptr<Rdm>(new BosRdm(exsig));
            else m_rdms[exsig] = std::unique_ptr<Rdm>(new FrmBosRdm(exsig));
        }

        for (auto exsig: m_specmom_exsigs) {
            REQUIRE_TRUE(exsig, "multidimensional estimators require a nonzero number of SQ operator indices");
            if (conn_utils::pure_frm(exsig))
                m_specmoms[exsig] = std::unique_ptr<SpectralMoment>(new FrmSpectralMoment(exsig, 1ul));
            else if (conn_utils::pure_bos(exsig))
                m_specmoms[exsig] = std::unique_ptr<SpectralMoment>(new BosSpectralMoment(exsig, 1ul));
            else m_specmoms[exsig] = std::unique_ptr<SpectralMoment>(new FrmBosSpectralMoment(exsig, 1ul));
        }
    }

    /**
     * should only be called if:
     *  the WF row is about to be removed
     *  the end of a sampling period has been reached
     *  the end of the calculation has been reached (all WF rows are about to be removed)
     * @param mbf
     */
    void make_diagonal_contribs(const fields::FrmOnv &mbf,
                                const size_t iroot_src, const size_t iroot_dst, const defs::wf_t& ci_product) {
        for (auto exsig: m_rdm_exsigs){
            m_rdms[exsig]->make_diagonal_contribs(mbf, iroot_src, iroot_dst, ci_product);
        }
    }


    defs::ham_comp_t estimate_energy(const Hamiltonian &ham) {
        return 0.0;
    }
};


#endif //M7_BILINEARESTIMATORGROUP_H
