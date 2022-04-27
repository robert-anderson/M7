//
// Created by Robert J. Anderson on 25/08/2021.
//

#ifndef M7_EXCITGENTESTER_H
#define M7_EXCITGENTESTER_H

#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/excititer/ExcitIter.h"

namespace excit_gen_tester {
    struct ResultRow : Row {
        field::MaeInds m_inds;
        field::Number<size_t> m_occur;
        field::Number<defs::prob_t> m_weight;
        field::Number<defs::ham_t> m_helem;

        ResultRow(size_t exsig) :
                m_inds(this, exsig, "excitation indices"),
                m_occur(this, "number of occurrences"),
                m_weight(this, "sum of reciprocal probabilities"),
                m_helem(this, "hamiltonian matrix element") {}

        field::MaeInds &key_field() {
            return m_inds;
        }
    };

    typedef BufferedTable<ResultRow, true> result_table_t;

    struct ExcitGenTester {
        ExcitGen &m_gen;
        ExcitIter &m_iter;
        result_table_t m_results;

        ExcitGenTester(ExcitGen &gen, ExcitIter &iter) :
                m_gen(gen), m_iter(iter), m_results("excit gen test results", {{iter.m_exsig}}) {
            m_results.set_expansion_factor(1.5);
            m_results.resize(1000);
        }

        /**
         * add a row to the results table for every outcome
         */
        template<typename mbf_t>
        void fill_results_table(const mbf_t &src_mbf) {
            typedef conn::from_field_t<mbf_t> conn_t;
            buffered::MaeInds work_inds(m_iter.m_exsig);
            m_results.clear();

            auto body_fn = [&](const conn_t &conn, defs::ham_t helem) {
                work_inds = conn;
                // if this key is already in the table then the iterator is emitting duplicate connections!
                DEBUG_ASSERT_EQ(*m_results[work_inds], ~0ul, "row should not already be mapped");
                auto irow = m_results.insert(work_inds);
                auto &row = m_results.m_row;
                row.jump(irow);
                row.m_occur = 0;
                row.m_weight = 0.0;
                row.m_helem = helem;
            };
            m_iter.foreach(src_mbf, body_fn, true);
        }

        template<typename mbf_t>
        size_t run(const mbf_t &src_mbf, size_t ndraw) {
            typedef conn::from_field_t<mbf_t> conn_t;
            CachedOrbs orbs(m_iter.m_ham.m_frm->m_point_group_map);
            conn_t conn(src_mbf);
            size_t exsig = m_iter.m_exsig;
            buffered::MaeInds work_inds(m_iter.m_exsig);
            size_t nnull = 0ul;

            DEBUG_ASSERT_FALSE(m_results.is_cleared(), "no connections were found by the excitation iterator");
            defs::prob_t prob = 0.0;
            defs::ham_t helem = 0.0;
            auto &row = m_results.m_row;
            for (size_t idraw = 0ul; idraw < ndraw; ++idraw) {
                auto success = m_gen.draw(exsig, src_mbf, orbs, prob, helem, conn);
                if (!success) {
                    ++nnull;
                    continue;
                }
                DEBUG_ASSERT_FALSE(consts::nearly_zero(prob), "non-null excitation generated with zero prob!");
                DEBUG_ASSERT_EQ(conn.exsig(), exsig, "generated excitation has the wrong exsig");
                work_inds = conn;
                auto irow = *m_results[work_inds];
                DEBUG_ASSERT_NE(irow, ~0ul, "excit generated that was not found in deterministic enumeration");
                row.jump(irow);
                row.m_occur++;
                row.m_weight += 1 / prob;
                DEBUG_ASSERT_TRUE(consts::nearly_equal(defs::ham_t(row.m_helem), helem),
                                  "excit gen returned the wrong H matrix element");
            }
            return nnull;
        }

        /**
         * @return
         *  true if the occurrence of each connection is non-zero
         */
        bool all_drawn_at_least_once() const;

        /**
         * check that all excitations were drawn with the correct probability
         * @param ndraw
         *  number of draws made
         * @param cutoff
         *  ratio of the occurrence to number of draws below which the normalized weight will not be considered
         * @param tol
         *  tolerance on which to decide correctness
         * @return
         *  true if all weights are correct within tolerance
         */
        bool all_correct_weights(size_t ndraw, double cutoff= 1e-2, defs::prob_t tol= 1e-2) const;

        /**
         * @param ndraw
         *  number of draws made
         * @return
         *  the average absolute error in the normalized weights
         */
        defs::prob_t mean_abs_error(size_t ndraw) const;
    };
}

#endif //M7_EXCITGENTESTER_H
