//
// Created by rja on 03/12/22.
//

#ifndef M7_COMMUTINGOBSERVABLE_H
#define M7_COMMUTINGOBSERVABLE_H

#include "M7_lib/wavefunction/Reference.h"

namespace commuting_obs {
        /**
     * commuting operators have simultaneous eigenfunctions, and so for operators A with the property [A, H] = 0
     * this can be exploited to obtain estimated eigenvalues from the (stochastic) wavefunction
     * A |Psi> = a |Psi>
     * i.e. < proj | A | Psi > / < proj | Psi > = a
     * where proj is a "trial wavefunction", which can either be a single reference MBF or a linear combination
     */
    struct Estimator {
        /**
         * operator commuting with H
         */
        const Hamiltonian* m_op;
        /**
         * pointer to the references
         */
        const wf::References* m_refs;

        NdReduction<ham_t, c_ndim_wf> m_proj_num;

        Estimator(const Hamiltonian* op, const wf::References* refs);

        /**
         * occupied MBFs connected to the reference must contribute to the numerator inner product <ref | H | mbf>
         * @param mbf
         * @param weights
         */
        void make_numerator_contribs(const Walker& walker);

        void begin_cycle(uint_t /*icycle*/);

        void end_cycle(uint_t /*icycle*/);

    };

    struct SpinSquare {
        const SpinSquareFrmHam m_frm_ham;
        const Hamiltonian m_op;
        Estimator m_est;
        SpinSquare(const sys::frm::Sector& sector, const wf::References* refs):
                m_frm_ham(sector), m_op(&m_frm_ham), m_est(&m_op, refs){}
    };
}


#endif //M7_COMMUTINGOBSERVABLE_H
