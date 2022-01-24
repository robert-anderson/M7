//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Dense.h"


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


dense::GemmWrapper::GemmWrapper(size_t nrowp, size_t ncolp, size_t nrowq, size_t ncolq,
                                char transp, char transq, size_t nrowr, size_t ncolr) :
        m_transa(valid_trans(transq)), m_transb(valid_trans(transp)),
        m_m(utils::safe_narrow<int>((m_transa=='N') ? ncolq : nrowq)),
        m_n(utils::safe_narrow<int>((m_transb=='N') ? nrowp : ncolp)),
        m_k(utils::safe_narrow<int>((m_transa=='N') ? nrowq : ncolq)),
        m_nrow_opb(utils::safe_narrow<int>((m_transb=='N') ? ncolp : nrowp)),
        m_lda(utils::safe_narrow<int>(ncolq)), m_ldb(utils::safe_narrow<int>(ncolp)),
        m_ldc(utils::safe_narrow<int>(ncolr)){
            REQUIRE_EQ(m_k, m_nrow_opb, "mismatch of contracted dimensions");
            REQUIRE_EQ(nrowr, (m_transb=='N') ? nrowp: ncolp, "incompatible result dimension");
            REQUIRE_EQ(ncolr, (m_transa=='N') ? ncolq: nrowq, "incompatible result dimension");
}

void dense::GemmWrapper::multiply(const float *p, const float *q, float *r, float alpha, float beta) {
    sgemm_(&m_transa, &m_transb, &m_m, &m_n, &m_k, &alpha, q, &m_lda, p, &m_ldb, &beta, r, &m_ldc);
}

void dense::GemmWrapper::multiply(const double *p, const double *q, double *r, double alpha, double beta) {
    dgemm_(&m_transa, &m_transb, &m_m, &m_n, &m_k, &alpha, q, &m_lda, p, &m_ldb, &beta, r, &m_ldc);
}

void dense::GemmWrapper::multiply(const std::complex<float> *p, const std::complex<float> *q, std::complex<float> *r,
                                  std::complex<float> alpha, std::complex<float> beta) {
    cgemm_(&m_transa, &m_transb, &m_m, &m_n, &m_k, &alpha, q, &m_lda, p, &m_ldb, &beta, r, &m_ldc);
}

void dense::GemmWrapper::multiply(const std::complex<double> *p, const std::complex<double> *q, std::complex<double> *r,
                                  std::complex<double> alpha, std::complex<double> beta) {
    zgemm_(&m_transa, &m_transb, &m_m, &m_n, &m_k, &alpha, q, &m_lda, p, &m_ldb, &beta, r, &m_ldc);
}

char dense::GemmWrapper::valid_trans(char t) {
    const std::string valid_chars = "NnCcTt";
    auto i = valid_chars.find(t);
    REQUIRE_LT(i, valid_chars.size(), "invalid character given");
    return valid_chars[2*(i/2)];
}