//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Dense.h"

template<>
void dense::Matrix<double>::mv_multiply(const std::vector<double> &in, std::vector<double> &out) {
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
void dense::Matrix<std::complex<double>>::mv_multiply(const std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out) {
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
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<double> work(lwork);
    int info;
    evals.resize(n);
    dsyev_(&jobz, &uplo, &n, mat.ptr(), &n, evals.data(), work.data(), &lwork, &info);
}

void dense::diag_inplace(SquareMatrix<std::complex<double>> &mat, std::vector<double> &evals) {
    const char jobz = 'V';
    const char uplo = 'U';
    const int n = mat.nrow();
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
    const int n = mat.nrow();
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
    for (size_t i=0ul; i<mat.nrow(); ++i) evals.push_back({wr[i], wi[i]});
}

dense::GemmWrapper::GemmWrapper(size_t nrowp, size_t ncolp, size_t nrowq, size_t ncolq) :
        m_m(utils::safe_narrow<int>(ncolq)),
        m_n(utils::safe_narrow<int>(nrowp)),
        m_k(utils::safe_narrow<int>(nrowq)),
        m_nrow_opb(utils::safe_narrow<int>(ncolp)){
    REQUIRE_EQ(m_k, m_nrow_opb, "mismatch of contracted dimensions");
}

dense::GemmWrapper::GemmWrapper(size_t nrowp, size_t ncolp, size_t nrowq, size_t ncolq, size_t nrowr, size_t ncolr) :
        GemmWrapper(nrowp, ncolp, nrowq, ncolq) {
    REQUIRE_EQ(nrowr, nrowp, "incompatible result dimension");
    REQUIRE_EQ(ncolr, ncolq, "incompatible result dimension");
}

void dense::GemmWrapper::multiply(const double *p, const double *q, double *r, double alpha, double beta) {
    dgemm_(&m_transa, &m_transb, &m_m, &m_n, &m_k, &alpha, q, &m_m, p, &m_nrow_opb, &beta, r, &m_m);
}

void dense::multiply(const dense::Matrix<double> &p, const dense::Matrix<double> &q,
                     dense::Matrix<double> &r, double alpha, double beta) {
    GemmWrapper wrapper(p.nrow(), p.ncol(), q.nrow(), q.ncol(), r.nrow(), r.ncol());
    wrapper.multiply(p.ptr(), q.ptr(), r.ptr(), alpha, beta);
}