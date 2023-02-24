//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Dense.h"

void dense::MatrixBase::set_sizes(uint_t nrow, uint_t ncol) {
    m_nrow = nrow;
    m_ncol = ncol;
    m_nelement = m_ncol * m_nrow;
    m_row_size = m_ncol * m_element_size;
    m_size = m_nrow * m_row_size;
}

bool dense::MatrixBase::compatible(const dense::MatrixBase &other) const {
    if (m_element_size != other.m_element_size) return false;
    return (m_nrow==other.m_nrow) && (m_ncol==other.m_ncol);
}

void dense::MatrixBase::resize(uint_t nrow, uint_t ncol) {
    auto no_copy_resize = [&](){
        set_sizes(nrow, ncol);
        m_bw.resize(m_nelement * m_element_size);
    };
    if (!m_nelement) {
        no_copy_resize();
        return;
    }
    const auto old = *this;
    no_copy_resize();
    nrow = std::min(m_nrow, old.m_nrow);
    ncol = std::min(m_ncol, old.m_ncol);

    if (i_can_globally_modify()) {
        for (uint_t irow = 0ul; irow < nrow; ++irow)
            std::copy(old.cbegin(irow), old.cbegin(irow) + ncol * m_element_size, begin(irow));
    }
    if (m_bw.node_shared()) mpi::barrier_on_node();
}

dense::MatrixBase::MatrixBase(uint_t nrow, uint_t ncol, uint_t element_size, bool node_shared) :
        m_buffer("", 1, node_shared), m_bw(&m_buffer), m_element_size(element_size) {
    resize(nrow, ncol);
}

dense::MatrixBase::MatrixBase(const dense::MatrixBase &other) : MatrixBase(other.m_nrow, other.m_ncol, other.m_element_size, other.m_bw.node_shared()){
    m_bw = other.m_bw;
}

dense::MatrixBase &dense::MatrixBase::operator=(const dense::MatrixBase &other) {
    DEBUG_ASSERT_TRUE(compatible(other), "cannot copy from incompatible matrix");
    resize(other.m_nrow, other.m_ncol);
    if (i_can_globally_modify()) m_bw = other.m_bw;
    return *this;
}

bool dense::MatrixBase::operator==(const dense::MatrixBase &other) const {
    return !std::memcmp(cbegin(), other.cbegin(), m_size);
}

void dense::MatrixBase::zero() {
    if (i_can_globally_modify()) m_bw.clear();
    m_bw.set_end(m_nrow);
}

void dense::MatrixBase::transpose() {
    auto tmp = *this;
    set_sizes(m_ncol, m_nrow);
    if (i_can_globally_modify()) {
        for (uint_t icol = 0ul; icol < m_ncol; ++icol) set_col(icol, tmp.cbegin(icol));
    }
    if (m_bw.node_shared()) mpi::barrier_on_node();
}

void dense::MatrixBase::reorder_rows(const uintv_t &order) {
    if (i_can_globally_modify())
        sort::reorder(begin(), m_row_size, order);
    if (m_bw.node_shared()) mpi::barrier_on_node();
}

void dense::MatrixBase::set(const void *src) {
    if (i_can_globally_modify()) std::memcpy(begin(), src, m_size);
    if (m_bw.node_shared()) mpi::barrier_on_node();
}

void dense::MatrixBase::set_row(uint_t irow, const void *src) {
    std::memcpy(begin(irow), src, m_row_size);
}

void dense::MatrixBase::get_row(uint_t irow, void *dst) const {
    std::memcpy(dst, cbegin(irow), m_row_size);
}

void dense::MatrixBase::set_col(uint_t icol, const void *src) {
    auto src_ptr = reinterpret_cast<const buf_t*>(src);
    for (uint_t irow = 0ul; irow < m_nrow; ++irow) {
        std::memcpy(begin(irow, icol), src_ptr, m_element_size);
        src_ptr += m_element_size;
    }
}

void dense::MatrixBase::get_col(uint_t icol, void *dst) {
    auto dst_ptr = reinterpret_cast<buf_t*>(dst);
    for (uint_t irow = 0ul; irow < m_nrow; ++irow) {
        std::memcpy(dst_ptr, cbegin(irow, icol), m_element_size);
        dst_ptr += m_element_size;
    }
}

dense::GemmWrapper::GemmWrapper(uint_t nrowp, uint_t ncolp, uint_t nrowq, uint_t ncolq,
                                char transp, char transq, uint_t nrowr, uint_t ncolr) :
        m_transa(valid_trans(transq)), m_transb(valid_trans(transp)),
        m_m(convert::safe_narrow<int>((m_transa=='N') ? ncolq : nrowq)),
        m_n(convert::safe_narrow<int>((m_transb=='N') ? nrowp : ncolp)),
        m_k(convert::safe_narrow<int>((m_transa=='N') ? nrowq : ncolq)),
        m_nrow_opb(convert::safe_narrow<int>((m_transb=='N') ? ncolp : nrowp)),
        m_lda(convert::safe_narrow<int>(ncolq)), m_ldb(convert::safe_narrow<int>(ncolp)),
        m_ldc(convert::safe_narrow<int>(ncolr)){
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
    const str_t valid_chars = "NnCcTt";
    auto i = valid_chars.find(t);
    REQUIRE_LT(i, valid_chars.size(), "invalid character given");
    return valid_chars[2*(i/2)];
}

bool dense::diag(const dense::SquareMatrix<float> &mat, v_t<float> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    v_t<float> work(lwork);
    int info;
    evals.resize(n);
    auto a = mat;
    ssyev_("N", "U", &n, a.tbegin(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<float> &mat, dense::SquareMatrix<float> &evecs, v_t<float> &evals) {
    REQUIRE_TRUE(mat.compatible(evecs), "non-compatibility between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    v_t<float> work(lwork);
    int info;
    evals.resize(n);
    evecs = mat;
    ssyev_("V", "U", &n, evecs.tbegin(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<float> &mat, v_t<std::complex<float>> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 4*n);
    v_t<float> work(lwork);
    int info;
    v_t<float> real_evals(n, 0.0);
    v_t<float> imag_evals(n, 0.0);
    auto a = mat;
    sgeev_("N", "N", &n, a.tbegin(), &n, real_evals.data(), imag_evals.data(),
           nullptr, &n, nullptr, &n, work.data(), &lwork, &info);
    arith::zip(real_evals, imag_evals, evals);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, v_t<double> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    v_t<double> work(lwork);
    int info;
    evals.resize(n);
    auto a = mat;
    dsyev_("N", "U", &n, a.tbegin(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, dense::SquareMatrix<double> &evecs, v_t<double> &evals) {
    REQUIRE_TRUE(mat.compatible(evecs), "non-compatibility between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 3 * n - 1);
    v_t<double> work(lwork);
    int info;
    evals.resize(n);
    evecs = mat;
    dsyev_("V", "U", &n, evecs.tbegin(), &n, evals.data(), work.data(), &lwork, &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<double> &mat, v_t<std::complex<double>> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 4*n);
    v_t<double> work(lwork);
    int info;
    v_t<double> real_evals(n, 0.0);
    v_t<double> imag_evals(n, 0.0);
    auto a = mat;
    dgeev_("N", "N", &n, a.tbegin(), &n, real_evals.data(), imag_evals.data(),
           nullptr, &n, nullptr, &n, work.data(), &lwork, &info);
    arith::zip(real_evals, imag_evals, evals);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<std::complex<double>> &mat, v_t<double> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 2 * n - 1);
    const int lrwork = std::max(1, 3 * n - 1);
    v_t<std::complex<double>> work(lwork);
    v_t<double> rwork(lrwork);
    int info;
    evals.resize(n);
    auto a = mat;
    zheev_("N", "U", &n, a.tbegin(), &n, evals.data(), work.data(), &lwork, rwork.data(), &info);
    return !info;
}

bool dense::diag(const dense::SquareMatrix<std::complex<double>> &mat, dense::SquareMatrix<std::complex<double>> &evecs,
                 v_t<double> &evals) {
    REQUIRE_TRUE(mat.compatible(evecs), "shape conflict between matrix and eigenvectors");
    const int n = mat.nrow();
    const int lwork = std::max(1, 2 * n - 1);
    const int lrwork = std::max(1, 3 * n - 1);
    v_t<std::complex<double>> work(lwork);
    v_t<double> rwork(lrwork);
    int info;
    evals.resize(n);
    evecs = mat;
    zheev_("V", "U", &n, evecs.tbegin(), &n, evals.data(), work.data(), &lwork, rwork.data(), &info);
    return !info;
}


bool dense::diag(const dense::SquareMatrix<std::complex<double>> &mat, v_t<std::complex<double>> &evals) {
    const int n = mat.nrow();
    const int lwork = std::max(1, 4*n);
    v_t<std::complex<double>> work(lwork);
    v_t<double> rwork(lwork);
    int info;
    evals.resize(n);
    auto a = mat;
    zgeev_("N", "N", &n, a.tbegin(), &n, evals.data(),
           nullptr, &n, nullptr, &n, work.data(), &lwork, rwork.data(), &info);
    return !info;
}
