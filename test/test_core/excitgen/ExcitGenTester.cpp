//
// Created by Robert J. Anderson on 25/08/2021.
//

#include "ExcitGenTester.h"

bool excit_gen_tester::ExcitGenTester::all_drawn_at_least_once() const {
    auto &row = m_results.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (row.m_occur == 0ul) return false;
    }
    return true;
}

bool excit_gen_tester::ExcitGenTester::all_correct_weights(uint_t ndraw, double cutoff, defs::prob_t tol) const {
    auto &row = m_results.m_row;
    for (row.restart(); row.in_range(); row.step()){
        auto ratio = double(row.m_occur)/ndraw;
        if (ratio<cutoff) continue;
        auto norm_weight = defs::prob_t(row.m_weight)/ndraw;
        if (std::abs(norm_weight-1.0) > tol) return false;
    }
    return true;
}

defs::prob_t excit_gen_tester::ExcitGenTester::mean_abs_error(uint_t ndraw) const {
    auto &row = m_results.m_row;
    defs::prob_t error = 0.0;
    for (row.restart(); row.in_range(); row.step()){
        auto norm_weight = defs::prob_t(row.m_weight)/ndraw;
        error += std::abs(norm_weight-1.0);
    }
    return error/m_results.m_hwm;
}
