//
// Created by rja on 19/02/23.
//

#ifndef M7_HFEXCITHISTS_H
#define M7_HFEXCITHISTS_H

#include "M7_lib/communication/SharedRows.h"
#include "M7_lib/mae/MaeTable.h"
#include "M7_lib/wavefunction/Wavefunction.h"

namespace hf_excit_hist {

    struct IndVals {
    private:
        uint_t m_nelement = 0ul;
    public:
        dense::Matrix<rdm_ind_t> m_inds;
        dense::Vector<wf_t> m_vals;
        IndVals(const hdf5::NodeReader& parent, str_t name, wf_t thresh);

        uint_t nelement() const {
            return m_nelement;
        }
    };

    struct Initializer {
        /**
         * wavefunction into which the permanitiators are to be inserted
         */
        wf::Fci& m_wf;
        /**
         * reference to Hartree--Fock (like) Many-body basis function
         */
        const field::Mbf& m_hf;
        /**
         * sorted indices and values of the C2 arrays
         */
        const IndVals m_c2;
        /**
         * minimum intermediate-normalized weight for inclusion as a permanitiator
         */
        const wf_comp_t m_thresh;
        /**
         * whether to revoke permanitiator status if multiple contributions sum to a value lower in magnitude than thresh
         */
        const bool m_cancellation;
        /**
         * working object in which C2^n |HF> is formed
         */
        buffered::Mbf m_work_mbf;
        /**
         * working object for computing the phase of excitations relative to the HF state
         */
        conn::Mbf m_work_conn;
        /**
         * number of permanitiators created by excitation level
         */
        reduction::NdArray<uint_t, 1> m_ncreated;

        Initializer(wf::Fci& wf, const field::Mbf& hf, str_t fname, wf_comp_t thresh, bool cancellation);

    private:

        /**
         * cumulatively multiply the first n ci coefficients until the product falls below the threshold
         */
        uint_t max_power();

        bool apply(field::Mbf& mbf, uint_t ientry);

        bool undo(field::Mbf& mbf, uint_t ientry);

        bool phase(field::Mbf& mbf);

        void setup(field::Mbf& mbf, uint_t imax, uint_t ipower, wf_t prev_product);

    public:
        void setup();
    };

    void initialize(wf::Fci& wf, const field::Mbf& hf, str_t fname, wf_t thresh, bool cancellation);

    void initialize(wf::Fci& wf, const field::Mbf& hf, const conf::CiPermanitiator& opts);

    struct Accumulators {
        /**
         * Hartree-Fock walker
         */
        const shared_rows::Walker* m_hf;
        /**
         * threshold for histogramming (intermediate normalization)
         */
        const wf_t m_thresh;
        /**
         * one table for each excitation level less than or equal to the user-specified max level
         */
        v_t<buffered::MappedTable<RdmRow>> m_tables;
        v_t<buffered::RdmInds> m_lookup_keys;
        /**
         * working object for computing connection to HF
         */
        mutable conn::Mbf m_work_conn;
        /**
         * average HF weight so that output can be intermediate normalized
         */
        wf_t m_norm = 0.0;
        const str_t m_save_file_name;

        const uintv_t m_nexcits;
        const uintv_t m_accumulated_nexcit_inds;

    private:
        static uintv_t make_nexcit_is_accumulated(const uintv_t& nexcits);

        uint_t ind(uint_t nexcit) const;

    public:
        Accumulators(const shared_rows::Walker* hf, uintv_t nexcits, wf_t thresh, str_t save_file_name);

        operator bool () const {
            return !m_tables.empty();
        }

        Accumulators(const conf::HfExcits& opts, const shared_rows::Walker* hf);

        ~Accumulators() {
            if (*this) save();
        }

        void add(const field::Mbf& mbf, wf_t weight);

        void save(const hdf5::NodeWriter& nw);

        void save();
    };
}

#endif //M7_HFEXCITHISTS_H
