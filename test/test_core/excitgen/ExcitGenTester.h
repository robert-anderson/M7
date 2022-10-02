//
// Created by Robert J. Anderson on 25/08/2021.
//

#ifndef M7_EXCITGENTESTER_H
#define M7_EXCITGENTESTER_H

#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/foreach/ConnForeach.h"
#include "M7_lib/hamiltonian/Hamiltonian.h"

namespace excit_gen_tester {
    struct ResultRow : Row {
        field::MaeInds m_inds;
        field::Number<uint_t> m_occur;
        field::Number<prob_t> m_weight;
        field::Number<ham_t> m_helem;

        ResultRow(uint_t exsig) :
                m_inds(this, exsig, "excitation indices"),
                m_occur(this, "number of occurrences"),
                m_weight(this, "sum of reciprocal probabilities"),
                m_helem(this, "hamiltonian matrix element") {}

        field::MaeInds &key_field() {
            return m_inds;
        }
    };

    typedef buffered::MappedTable<ResultRow> result_table_t;


    struct ExcitGenTester {
        const Hamiltonian &m_h;
        ExcitGen& m_excit_gen;
        conn_foreach::Base &m_conn_iter;
        result_table_t m_results;

        /**
         * a success flag, and flags corresponding to each of the reasons for which an excitation generator can fail
         */
        enum Status {
            Success, /* no issues with draw loop */
            AllNull, GenWithZeroProb, InvalidDrawProb, InvalidProb, InvalidProbGivenHElem,
            ProbMismatch, ProbMismatchGivenHElem, WrongExsig, Unconnected, WrongHElem,
            NotAllDrawnAtLeastOnce, VarianceIncrease, WrongWeights
        };

        /**
         * each non-success Status is associated with an error message
         * @param status
         *  status flag emitted by the draw loop
         * @return
         *  associated error message
         */
        str_t status_error_msg(Status status);

        ExcitGenTester(const Hamiltonian &h, ExcitGen &excit_gen, conn_foreach::Base &conn_iter) :
                m_h(h), m_excit_gen(excit_gen), m_conn_iter(conn_iter),
                m_results("excit gen test results", {{conn_iter.m_exsig}}) {
            m_results.set_expansion_factor(1.5);
            m_results.resize(1000);
        }

        /**
         * add a row to the results table for every outcome
         */
        template<typename mbf_t>
        void fill_results_table(const mbf_t &src_mbf) {
            typedef conn::from_field_t<mbf_t> conn_t;
            buffered::MaeInds work_inds(m_conn_iter.m_exsig);
            m_results.clear();
            conn_t conn(m_h.m_basis.size());
            auto body_fn = [&]() {
                auto helem = m_h.get_element(src_mbf, conn);
                if (!ham::is_significant(helem)) return;
                work_inds = conn;
                // if this key is already in the table then the iterator is emitting duplicate connections!
                DEBUG_ASSERT_TRUE(m_results.lookup(work_inds), "row should not already be mapped");
                auto irow = m_results.insert(work_inds);
                auto &row = m_results.m_row;
                row.jump(irow);
                row.m_occur = 0;
                row.m_weight = 0.0;
                row.m_helem = helem;
            };
            m_conn_iter.loop(conn, src_mbf, body_fn);
        }

        template<typename mbf_t>
        Status perform_draws(const mbf_t &src_mbf, uint_t ndraw) {
            typedef conn::from_field_t<mbf_t> conn_t;
            conn_t conn(src_mbf);
            uint_t exsig = m_conn_iter.m_exsig;
            buffered::MaeInds work_inds(m_conn_iter.m_exsig);
            uint_t nnull = 0ul;

            REQUIRE_FALSE(m_results.is_cleared(), "no connections were found by the excitation iterator");
            prob_t prob = 0.0;
            ham_t helem = 0.0;
            auto &row = m_results.m_row;
            for (uint_t idraw = 0ul; idraw < ndraw; ++idraw) {

                auto success = m_excit_gen.draw(exsig, src_mbf, prob, helem, conn);
                if (!success) {
                    ++nnull;
                    continue;
                }

                if (prob < 0.0 || prob > 1.0) {
                    return InvalidDrawProb;
                }
                if (fptol::numeric_zero(prob)) return GenWithZeroProb;
                const auto chk_prob = m_excit_gen.prob(src_mbf, conn);
                if (chk_prob < 0.0 || chk_prob > 1.0) return InvalidProb;
                if (!fptol::numeric_equal(prob, chk_prob)) {
                    return ProbMismatch;
                }
                const auto chk_prob_given_helem = m_excit_gen.prob(src_mbf, conn, helem);
                if (chk_prob_given_helem < 0.0 || chk_prob_given_helem > 1.0) return InvalidProbGivenHElem;
                if (!fptol::numeric_equal(prob, chk_prob_given_helem)) return ProbMismatchGivenHElem;
                if (conn.exsig()!=exsig) return WrongExsig;
                work_inds = conn;
                auto lookup = m_results.lookup(work_inds, row);
                if (!lookup) return Unconnected;
                row.m_occur++;
                row.m_weight += 1 / prob;
                if (!fptol::numeric_equal(ham_t(row.m_helem), helem)) return WrongHElem;
            }
            return Success;
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
         *  true if all weights are correct within absolute tolerance
         */
        bool all_correct_weights(uint_t ndraw, double cutoff, prob_t atol) const;

        /**
         * @param ndraw
         *  number of draws made
         * @return
         *  the average absolute error in the normalized weights
         */
        prob_t mean_abs_error(uint_t ndraw) const;

    private:
        template<typename mbf_t>
        void run(const mbf_t &src_mbf, uint_t ndraw, Status& status, double cutoff, prob_t atol) {
            fill_results_table(src_mbf);
            status = perform_draws(src_mbf, ndraw);
            if (status!=Success) return;
            status = all_drawn_at_least_once() ? Success : NotAllDrawnAtLeastOnce;
            if (status!=Success) return;
            auto av_err1 = mean_abs_error(ndraw);
            status = perform_draws(src_mbf, ndraw);
            if (status!=Success) return;
            auto av_err2 = mean_abs_error(2 * ndraw);
            status = (av_err2 < av_err1) ? Success : VarianceIncrease;
            if (status!=Success) return;
            status = all_correct_weights(2 * ndraw, cutoff, atol) ? Success : WrongWeights;
        }
    public:
        template<typename mbf_t>
        std::string run(const mbf_t &src_mbf, uint_t ndraw, double cutoff = 1e-2, prob_t atol = 1e-2) {
            Status status;
            run(src_mbf, ndraw, status, cutoff, atol);
            return status_error_msg(status);
        }
    };
}

#endif //M7_EXCITGENTESTER_H
