//
// Created by rja on 03/03/2020.
//

#include "EigenSolver.h"

template<>
void EigenSolver<double>::execute() {
    dsyev_(&m_jobz, &m_uplo, &m_n, &m_evecs(0, 0), &m_n, m_evals.data(), m_work.data(), &m_lwork, &m_info);
}

template<>
void EigenSolver<std::complex<double>>
::execute() {
    if (!m_rwork)
        m_rwork = std::unique_ptr<std::vector<double>>(new std::vector<double>(std::max(1, 3 * m_n - 2)));
    zheev_(&m_jobz, &m_uplo, &m_n, &m_evecs(0, 0), &m_n, m_evals.data(), m_work.data(), &m_lwork, m_rwork->data(), &m_info
    );
}
