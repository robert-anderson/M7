//
// Created by rja on 19/02/23.
//

#ifndef M7_HFEXCITHISTS_H
#define M7_HFEXCITHISTS_H

#include "M7_lib/communication/SharedRows.h"
#include "M7_lib/mae/MaeTable.h"

namespace hf_excit_hist {

    struct Initializer {

    };

    struct Initializers {

    };

//    struct Accumulator {
//        buffered::MappedTable<RdmRow> m_table;
//        void save(const hdf5::NodeWriter& nw) const {
//        }
//        void save() {
//        }
//    };

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
        wf_t m_norm;
        const str_t m_save_file_name;

        Accumulators(const shared_rows::Walker* hf, uint_t nexcit_max, wf_t thresh, str_t save_file_name):
            m_hf(hf), m_thresh(thresh), m_work_conn(hf->mbf()), m_save_file_name(std::move(save_file_name)){
            if (!nexcit_max) return;
            REQUIRE_TRUE(m_hf, "HF state must be defined for HF excitation accumulation");
            m_tables.reserve(nexcit_max);
            m_lookup_keys.reserve(nexcit_max);
            for (uint_t nexcit = 1ul; nexcit < nexcit_max + 1; ++nexcit) {
                const OpSig exsig({nexcit, nexcit}, {0, 0});
                const auto name = exsig.to_string()+" excitations of HF state";
                m_tables.emplace_back(name, RdmRow(exsig, 1), false);
                m_lookup_keys.emplace_back(exsig);
            }
        }

        operator bool () const {
            return !m_tables.empty();
        }

        Accumulators(const conf::HfExcits& opts, const shared_rows::Walker* hf):
            Accumulators(hf, opts.m_max_nexcit, opts.m_thresh, opts.m_save.m_path){}

        ~Accumulators() {
            if (*this) save();
        }

        void add(const field::Mbf& mbf, wf_t weight) {
            if (!*this) return;
            m_work_conn.connect(m_hf->mbf(), mbf);
            const auto nexcit = m_work_conn.exsig().nfrm_cre();
            if (nexcit > m_tables.size() + 1) return;
            if (!nexcit) {
                m_norm += weight;
                return;
            }
            auto& table = m_tables[nexcit - 1];
            table.m_row.m_inds = m_work_conn;
            auto& lookup = table.lookup(table.m_row.m_inds);
            // if the excitation is already in the table, it is added regardless of current weight
            if (lookup) lookup.m_values[0] += weight;
            // if it's not already histogrammed, it can be added only if the instantaneous weight is sufficient
            else if (weight >= m_hf->weight(0)*m_thresh) {
                table.m_row.m_values = weight;
                table.insert(table.m_row);
            }
        }

        void save(const hdf5::NodeWriter& nw) {
            uint_t i=1ul;
            for (auto& table : m_tables) {
                auto& row = table.m_row;
                for (row.restart(); row; ++row) row.m_values /= m_norm;
                table.save(nw, OpSig({i, i}, {0, 0}).to_string(), true);
                ++i;
            }
        }

        void save() {
            hdf5::FileWriter fw(m_save_file_name);
            save(fw);
        }
    };
}

#endif //M7_HFEXCITHISTS_H
