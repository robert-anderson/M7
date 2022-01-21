//
// Created by rja on 03/03/2020.
//

#ifndef M7_EIGENSOLVER_H
#define M7_EIGENSOLVER_H

#include "Dense.h"

template<typename T>
class EigenSolver {
    char m_jobz = 'V';
    char m_uplo = 'U';
    int m_n;
    int m_lwork;
    int m_info;
    std::vector<T> m_work;
    std::unique_ptr<std::vector<consts::comp_t<T>>> m_rwork = nullptr;
public:
    dense::Matrix<T> m_evecs;
    std::vector<consts::comp_t<T>> m_evals;

    EigenSolver(const dense::Matrix<T> &matrix) :
            m_n((int) matrix.nrow()), m_lwork(std::max(1, 3 * m_n - 1)), m_work(m_lwork),
            m_evecs(matrix), m_evals(matrix.nrow()) {
        execute();
    }

private:
    void execute();
};


#endif //M7_EIGENSOLVER_H
