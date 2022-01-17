//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Dense.h"

template<>
void dense::Matrix<double>::multiply(const std::vector<double> &in, std::vector<double> &out) {
    char trans = 'N';
    const double alpha = 1.0;
    const double beta = 1.0;
    int nrow = m_nrow;
    int ncol = m_ncol;
    int incx = 1;
    int incy = 1;
    dgemv_(&trans, &ncol, &nrow, &alpha, m_buffer.data(), &nrow, in.data(), &incx, &beta, out.data(), &incy);
}

template<>
void dense::Matrix<std::complex<double>>::multiply(const std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out) {
    char trans = 'N';
    const std::complex<double> alpha = 1.0;
    const std::complex<double> beta = 1.0;
    int nrow = m_nrow;
    int ncol = m_ncol;
    int incx = 1;
    int incy = 1;
    zgemv_(&trans, &ncol, &nrow, &alpha, m_buffer.data(), &nrow, in.data(), &incx, &beta, out.data(), &incy);
}

void dense::diag_inplace(SquareMatrix<double> &mat, std::vector<double> &evals) {
    const char jobz = 'V';
    const char uplo = 'U';
    const int n = mat.m_nrow;
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<double> work(lwork);
    int info;
    evals.resize(n);
    dsyev_(&jobz, &uplo, &n, mat.ptr(), &n, evals.data(), work.data(), &lwork, &info);
}

void dense::diag_inplace(SquareMatrix<std::complex<double>> &mat, std::vector<double> &evals) {
    const char jobz = 'V';
    const char uplo = 'U';
    const int n = mat.m_nrow;
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(std::max(1, 3 * n - 2));
    int info;
    evals.resize(n);
    zheev_(&jobz, &uplo, &n, mat.ptr(), &n, evals.data(), work.data(), &lwork, rwork.data(), &info);
}

void dense::diagonalize(SquareMatrix<double>& mat, std::vector<std::complex<double>>& evals) {
    char jobvl = 'N';
    char jobvr = 'N';
    const int n = mat.m_nrow;
    std::vector<double> wr(n);
    std::vector<double> wi(n);
    const int ldvl = 1;
    const int ldvr = 1;
    const int lwork = std::max(1, 3 * n);
    int info;
    std::vector<double> work(lwork);
    dgeev_(&jobvl, &jobvr, &n, mat.ptr(), &n, wr.data(), wi.data(), nullptr,
           &ldvl, nullptr, &ldvr, work.data(), &lwork, &info);
    evals.clear();
    evals.reserve(n);
    for (size_t i=0ul; i<mat.m_nrow; ++i) evals.push_back({wr[i], wi[i]});
}
