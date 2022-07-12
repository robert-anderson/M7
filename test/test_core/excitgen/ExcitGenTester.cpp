//
// Created by Robert J. Anderson on 25/08/2021.
//

#include "ExcitGenTester.h"


str_t excit_gen_tester::ExcitGenTester::status_error_msg(Status status) {
    switch (status) {
        case Success: return "";
        case AllNull: return "no non-null excitations were drawn";
        case GenWithZeroProb: return "non-null excitation generated with zero probability";
        case InvalidDrawProb: return "draw probability is not in valid range";
        case InvalidProb: return "probability is not in valid range";
        case InvalidProbGivenHElem:
            return "probability given the H matrix element is not in valid range";
        case ProbMismatch: return "prob of connection doesn't match prob resulting from the draw method";
        case ProbMismatchGivenHElem:
            return "prob of connection doesn't match prob resulting from the draw method given H matrix element";
        case WrongExsig: return "drawn connection has wrong excitation signature";
        case Unconnected: return "excit generated that was not found in deterministic enumeration";
        case WrongHElem: return "H element returned by draw method is not consistent with value from Hamiltonian class";
        case NotAllDrawnAtLeastOnce: return "not all expected connections were drawn at least once";
        case VarianceIncrease: return "variance in normalized weights did not decrease after doubling number of draws";
        case WrongWeights: return "normalized weights are not all correct within given tolerance";
    }
    return "unknown error";
}

bool excit_gen_tester::ExcitGenTester::all_drawn_at_least_once() const {
    auto &row = m_results.m_row;
    for (row.restart(); row.in_range(); row.step()){
        if (row.m_occur == 0ul) return false;
    }
    return true;
}

bool excit_gen_tester::ExcitGenTester::all_correct_weights(uint_t ndraw, double cutoff, prob_t tol) const {
    auto &row = m_results.m_row;
    for (row.restart(); row.in_range(); row.step()){
        auto ratio = double(row.m_occur)/ndraw;
        if (ratio<cutoff) continue;
        auto norm_weight = prob_t(row.m_weight)/ndraw;
        if (std::abs(norm_weight-1.0) > tol) return false;
    }
    return true;
}

prob_t excit_gen_tester::ExcitGenTester::mean_abs_error(uint_t ndraw) const {
    auto &row = m_results.m_row;
    prob_t error = 0.0;
    for (row.restart(); row.in_range(); row.step()){
        auto norm_weight = prob_t(row.m_weight)/ndraw;
        error += std::abs(norm_weight-1.0);
    }
    return error/m_results.m_hwm;
}