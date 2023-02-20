//
// Created by rja on 19/02/23.
//

#ifndef M7_HFEXCITHISTS_H
#define M7_HFEXCITHISTS_H

#include "M7_lib/communication/SharedRows.h"
#include "M7_lib/mae/MaeTable.h"

namespace hf_excit_hist {

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

        const uintv_t m_nexcits;
        const uintv_t m_accumulated_nexcit_inds;

    private:
        static uintv_t make_nexcit_is_accumulated(const uintv_t& nexcits) {
            if (nexcits.empty()) return {};
            auto nmax = *std::max_element(nexcits.cbegin(), nexcits.cend());
            uintv_t out(nmax+1, ~0ul);
            uint_t i = 0;
            for (auto& n: nexcits) out[n] = i++;
            return out;
        }

        uint_t ind(uint_t nexcit) const {
            return nexcit < m_accumulated_nexcit_inds.size() ? m_accumulated_nexcit_inds[nexcit] : ~0ul;
        }

    public:
        Accumulators(const shared_rows::Walker* hf, uintv_t nexcits, wf_t thresh, str_t save_file_name):
            m_hf(hf), m_thresh(thresh), m_work_conn(hf->mbf()), m_save_file_name(std::move(save_file_name)),
            m_nexcits(std::move(nexcits)), m_accumulated_nexcit_inds(make_nexcit_is_accumulated(m_nexcits)) {
            if (m_nexcits.empty()) return;
            REQUIRE_TRUE(m_hf, "HF state must be defined for HF excitation accumulation");
            m_tables.reserve(m_nexcits.size());
            m_lookup_keys.reserve(m_nexcits.size());
            for (auto& nexcit: m_nexcits){
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
            Accumulators(hf, opts.m_nexcits, opts.m_thresh, opts.m_save.m_path){}

        ~Accumulators() {
            if (*this) save();
        }

        void add(const field::Mbf& mbf, wf_t weight) {
            if (!*this) return;
            m_work_conn.connect(m_hf->mbf(), mbf);
            const auto nexcit = m_work_conn.exsig().nfrm_cre();
            if (!nexcit) {
                m_norm += weight;
                return;
            }
            const auto i = ind(nexcit);
            if (i == ~0ul) return;
            auto& table = m_tables[i];
            auto& key = m_lookup_keys[i];
            key = m_work_conn;
            auto& lookup = table.lookup(key);
            // if the excitation is already in the table, it is added regardless of current weight
            if (lookup) lookup.m_values[0] += weight;
            // if it's not already histogrammed, it can be added only if the instantaneous weight is sufficient
            else if (std::abs(weight) >= std::abs(m_hf->weight(0)*m_thresh)) {
                auto& insert_row = table.insert(key);
                insert_row.m_values += weight;
            }
        }

        void save(const hdf5::NodeWriter& nw) {
            auto table = m_tables.begin();
            for (auto nexcit: m_nexcits) {
                OpSig exsig({nexcit, nexcit}, {0, 0});
                auto& row = table->m_row;
                for (row.restart(); row; ++row) row.m_values /= m_norm;
                table->save(nw, exsig.to_string(), true);
                ++table;
            }
        }

        void save() {
            hdf5::FileWriter fw(m_save_file_name);
            save(fw);
        }
    };
}

#endif //M7_HFEXCITHISTS_H
