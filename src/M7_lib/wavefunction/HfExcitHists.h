//
// Created by rja on 19/02/23.
//

#ifndef M7_HFEXCITHISTS_H
#define M7_HFEXCITHISTS_H

#include "M7_lib/communication/SharedRows.h"
#include "M7_lib/mae/MaeTable.h"
#include "M7_lib/parallel/PeriodicEvent.h"

namespace wf {
    struct Vectors;
}

namespace hf_excit_hist {
    struct IndVals {
        wf_comp_t m_geo_mean = 0.0;
        dense::Matrix<rdm_ind_t> m_inds;
        dense::Vector<wf_t> m_vals;

        IndVals(const hdf5::NodeReader& parent, str_t name);

        uint_t nelement() const {
            return m_vals.nelement();
        }
    };

    struct Initializer {
        /**
         * wavefunction into which the permanitiators are to be inserted
         */
        wf::Vectors& m_wf;
        /**
         * reference to Hartree--Fock (like) Many-body basis function
         */
        const field::Mbf& m_hf;
        /**
         * sorted indices and values of the C2 arrays
         */
        const IndVals m_c2;
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
        /**
         * k_i for each power i
         */
        const v_t<double> m_min_ks;
        /**
         * k_i + delta k_i for each power i
         */
        const v_t<double> m_max_ks;
        /**
         * C2 geometric mean to the power of max_ks for each power i
         */
        const v_t<double> m_threshs;

        Initializer(wf::Vectors& wf, const field::Mbf& hf, str_t fname,
                    uint_t max_exlvl, double delta_k, bool cancellation);

    private:
        /**
         * @return
         *  threshold which would yield a single permanitiator for C2^ipower excitations
         */
        wf_comp_t thresh_for_first_pmntr(uint_t ipower);

        /**
         * @return
         *  power of the C2 geometric mean corresponding to thresh_for_first_pmntr(ipower)
         */
        wf_comp_t gmp_for_first_pmntr(uint_t ipower);

        v_t<double> make_min_ks(uint_t npower);

        v_t<double> make_max_ks(double delta_k);

        v_t<double> make_threshs();

        bool apply(field::FrmOnv& mbf, uint_t ientry);
        bool apply(field::BosOnv&, uint_t) {return false;}
        bool apply(field::FrmBosOnv& mbf, uint_t ientry) {return apply(mbf.m_frm, ientry);}

        bool undo(field::FrmOnv& mbf, uint_t ientry);
        bool undo(field::BosOnv&, uint_t) {return false;}
        bool undo(field::FrmBosOnv& mbf, uint_t ientry) {return undo(mbf.m_frm, ientry);}

        template<typename mbf_t>
        bool phase(mbf_t& mbf) {
            m_work_conn.connect(m_hf, mbf);
            return m_work_conn.phase(m_hf);
        }

        /**
         * recurse through coefficient tree from ielement
         */
        bool loop_body(field::Mbf& mbf, uint_t ielement, uint_t ipower, wf_t prev_product);

        /**
         * generic loop
         */
        void setup(field::Mbf& mbf, uint_t imax, uint_t ipower, wf_t prev_product);

        /**
         * top-level loop
         */
        void setup(field::Mbf& mbf);

        void communicate_and_insert();

    public:
        /**
         * calls top-level loop and applies cancellation if required
         */
        void setup();
    };

    void initialize(wf::Vectors& wf, const field::Mbf& hf, str_t fname,
                    uint_t max_exlvl, double delta_k, bool cancellation);

    void initialize(wf::Vectors& wf, const field::Mbf& hf, const conf::CiPmntr& opts);

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
         * average HF weight
         */
        wf_t m_norm = 0.0;

        const uintv_t m_nexcits;
        const uintv_t m_accumulated_nexcit_inds;

        Epoch m_accum_epoch;
        const str_t m_save_file_path;
        PeriodicFileSeries m_chkpt_files;


    private:
        static uintv_t make_nexcit_is_accumulated(const uintv_t& nexcits);

        uint_t ind(uint_t nexcit) const;

    public:
        Accumulators(const shared_rows::Walker* hf, uintv_t nexcits, wf_t thresh,
                     const conf::OptionalFile& save_file,
                     const conf::OptionalFileSeries& chkpt_files);

        operator bool () const {
            return !m_tables.empty();
        }

        Accumulators(const conf::HfExcits& opts, const shared_rows::Walker* hf);

        ~Accumulators() {
            if (*this) save();
        }

        void add(const field::Mbf& mbf, wf_t weight);

        void save(const hdf5::NodeWriter& nw) const;

        void save() const;

        void attempt_chkpt(uint_t icycle);
    };
}

#endif //M7_HFEXCITHISTS_H
