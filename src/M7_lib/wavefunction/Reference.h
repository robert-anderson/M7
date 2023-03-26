//
// Created by Robert J. Anderson on 03/07/2020.
//

#ifndef M7_REFERENCE_H
#define M7_REFERENCE_H

#include <M7_lib/observables/HfExcits.h>
#include <M7_lib/basis/Suites.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/parallel/Reduction.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/nd/NdArray.h>
#include <M7_lib/communication/SharedRows.h>

#include "WalkerTable.h"

namespace wf {

    struct Vectors;

    class Ref : public shared_rows::Walker {
        using shared_rows::Walker::row;
        Vectors& m_wf;
        /**
         * index to the "part" of the wavefunction for which this object tracks the reference row
         */
        const uint_t m_ipart;
        /**
         * If a candidate for redefinition of the reference is found, then its weight and row within m_list must be stored
         */
        uint_t m_irow_candidate;
        /**
         * weight of current best candidate on which to redefine the reference row
         */
        wf_comp_t m_candidate_weight = 0.0;
        /**
         * default scale factor defining when the candidate weight has grown to the magnitude at which it must be made the
         * new reference for this WF part
         */
        const double m_redefinition_thresh;

        reduction::Scalar<ham_t> m_proj_energy_num;

    public:
        Ref(const conf::Reference& opts, Vectors& wf, uint_t ipart, TableBase::Loc loc);

        uint_t occupied_ncycle(uint_t icycle) const {
            return gathered().m_row.occupied_ncycle(icycle);
        }

        const wf_t& weight() const {
            return Walker::weight(m_ipart);
        }

        ham_comp_t hdiag() const {
            return gathered().m_row.m_hdiag;
        }

        /**
         * If the candidate with the highest absolute weight across all MPI ranks has a weight larger than that of the current
         * reference by a specified factor, then accept the candidate and let both the dynamic row set and reference
         * connection flags in the WF m_store table reflect this change
         * @param redefinition_thresh
         *  scale factor by which the candidate weight must exceed that of the current ref
         */
        void accept_candidate(uint_t icycle);

        /**
         * add contributions from the current m_wf.m_store.m_row
         */
        void contrib_row(const ::Walker& walker);

        /**
         * reset variables to begin a fresh MC cycle
         */
        void begin_cycle(uint_t icycle);

        void end_cycle(uint_t /*icycle*/);

        /**
         * @param mbf
         *  reference to MBF
         * @return
         *  true if the matrix element of the Hamiltonian between the current reference and the argument is non-zero
         */
        bool connected(const Mbf& mbf) const;

        /**
         * occupied MBFs connected to the reference must contribute to the numerator inner product <ref | H | mbf>
         * @param mbf
         * @param weights
         */
        void make_numerator_contribs(const Mbf& mbf, const wf_t& weight);

        const ham_t& proj_energy_num() const;

    };

    /**
     * Wavefunctions have many parts in general, and each of these may require different reference ONVs, so here we
     * define a vector of the above-defined single part Ref class
     */
    struct Refs {
        v_t<Ref> m_refs;
        buffered::Numbers<ham_t, c_ndim_wf> m_proj_energy_nums;
        buffered::Numbers<wf_t, c_ndim_wf> m_weights;

        Refs(const conf::Reference& opts, Vectors& wf, v_t<TableBase::Loc> locs);

        const Ref& operator[](const uint_t& ipart) const;

        void begin_cycle(uint_t icycle);

        void end_cycle(uint_t icycle);

        void contrib_row(const Walker& walker);

        v_t<bool> connected(const Mbf& mbf) const;

        const Numbers<ham_t, c_ndim_wf>& proj_energy_nums();

        const Numbers<wf_t, c_ndim_wf>& weights();

    };
}

#endif //M7_REFERENCE_H
