//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Dense.h"


dense::GemmWrapper::GemmWrapper(size_t nrowp, size_t ncolp, size_t nrowq, size_t ncolq,
                                char transp, char transq, size_t nrowr, size_t ncolr) :
        m_transa(valid_trans(transq)), m_transb(valid_trans(transp)),
        m_m(utils::convert::safe_narrow<int>((m_transa=='N') ? ncolq : nrowq)),
        m_n(utils::convert::safe_narrow<int>((m_transb=='N') ? nrowp : ncolp)),
        m_k(utils::convert::safe_narrow<int>((m_transa=='N') ? nrowq : ncolq)),
        m_nrow_opb(utils::convert::safe_narrow<int>((m_transb=='N') ? ncolp : nrowp)),
        m_lda(utils::convert::safe_narrow<int>(ncolq)), m_ldb(utils::convert::safe_narrow<int>(ncolp)),
        m_ldc(utils::convert::safe_narrow<int>(ncolr)){
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

bool dense::diag(const dense::SquareMatrix<float> &mat, std::vector<float> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<float> work(lwork);
    int info;
    evals.resize(n);
    auto a = mat;
    ssyev_("N", "U", &n, a.ptr(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<float> &mat, dense::SquareMatrix<float> &evecs, std::vector<float> &evals) {
    REQUIRE_TRUE(mat.dims()==evecs.dims(), "shape conflict between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<float> work(lwork);
    int info;
    evals.resize(n);
    evecs = mat;
    ssyev_("V", "U", &n, evecs.ptr(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<float> &mat, std::vector<std::complex<float>> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 4*n);
    std::vector<float> work(lwork);
    int info;
    std::vector<float> real_evals(n, 0.0);
    std::vector<float> imag_evals(n, 0.0);
    auto a = mat;
    sgeev_("N", "N", &n, a.ptr(), &n, real_evals.data(), imag_evals.data(),
           nullptr, &n, nullptr, &n, work.data(), &lwork, &info);
    complex_utils::combine(real_evals, imag_evals, evals);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, std::vector<double> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<double> work(lwork);
    int info;
    evals.resize(n);
    auto a = mat;
    dsyev_("N", "U", &n, a.ptr(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, dense::SquareMatrix<double> &evecs, std::vector<double> &evals) {
    REQUIRE_TRUE(mat.dims()==evecs.dims(), "shape conflict between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    std::vector<double> work(lwork);
    int info;
    evals.resize(n);
    evecs = mat;
    dsyev_("V", "U", &n, evecs.ptr(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, std::vector<std::complex<double>> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 4*n);
    std::vector<double> work(lwork);
    int info;
    std::vector<double> real_evals(n, 0.0);
    std::vector<double> imag_evals(n, 0.0);
    auto a = mat;
    dgeev_("N", "N", &n, a.ptr(), &n, real_evals.data(), imag_evals.data(),
           nullptr, &n, nullptr, &n, work.data(), &lwork, &info);
    complex_utils::combine(real_evals, imag_evals, evals);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<std::complex<double>> &mat, std::vector<double> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 2 * n - 1);
    const int lrwork = std::max(1, 3 * n - 1);
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(lrwork);
    int info;
    evals.resize(n);
    auto a = mat;
    zheev_("N", "U", &n, a.ptr(), &n, evals.data(), work.data(), &lwork, rwork.data(), &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<std::complex<double>> &mat, dense::SquareMatrix<std::complex<double>> &evecs,
                 std::vector<double> &evals) {
    REQUIRE_TRUE(mat.dims()==evecs.dims(), "shape conflict between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 2 * n - 1);
    const int lrwork = std::max(1, 3 * n - 1);
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> rwork(lrwork);
    int info;
    evals.resize(n);
    evecs = mat;
    zheev_("V", "U", &n, evecs.ptr(), &n, evals.data(), work.data(), &lwork, rwork.data(), &info);
    return !info;
}
