//
// Created by Robert John Anderson on 2020-01-24.
//

#ifndef M7_DENSE_H
#define M7_DENSE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cstring>
#include <src/core/parallel/MPIAssert.h>
#include "src/core/util/consts.h"
#include "src/defs.h"
#include "Sparse.h"


template<typename T>
class EigenSolver;

extern "C" void dsyev_(const char *jobz, const char *uplo, const int *n, double *a, const int *lda, double *w,
                       double *work, const int *lwork, int *info);

extern "C" void zheev_(const char *jobz, const char *uplo, const int *n, std::complex<double> *a, const int *lda,
                       double *w, std::complex<double> *work, const int *lwork, double *rwork, int *info);


extern "C" void dgeev_(const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda,
                       double *wr, double *wi,
                       double *vl, const int* ldvl,
                       double *vr, const int* ldvr,
                       double *work, const int *lwork, int *info);


extern "C" void dgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const double *alpha,
                       const double *a, const int *lda,
                       const double *b, const int *ldb,
                       const double *beta,
                       double *c, const int *ldc);

extern "C" void dgemv_(const char *trans, const int *m, const int *n,
                       const double *alpha, const double *a, const int *lda,
                       const double *x, const int *incx, const double *beta,
                       double *y, const int *incy);

extern "C" void zgemv_(const char *trans, const int *m, const int *n, const std::complex<double> *alpha,
                       const std::complex<double> *a, const int *lda, const std::complex<double> *x, const int *incx,
                       const std::complex<double> *beta, std::complex<double> *y, const int *incy);

namespace dense {

    template<typename T>
    class Matrix {
        std::vector<T> m_buffer;

        size_t index(const size_t &irow, const size_t &icol) const {
            DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
            DEBUG_ASSERT_LT(icol, m_ncol, "column index OOB");
            return irow * m_ncol + icol;
        }

        size_t m_nrow, m_ncol;
    public:

        Matrix(size_t nrow, size_t ncol) : m_buffer(nrow * ncol, T(0)), m_nrow(nrow), m_ncol(ncol) {}
        Matrix(const sparse::Matrix<T>& sparse) : Matrix(sparse.nrow(), sparse.max_column_index()+1){
            *this = sparse;
        }

        Matrix(const Matrix& other): Matrix(other.m_nrow, other.m_ncol){
            m_buffer = other.m_buffer;
        }
        Matrix& operator=(const Matrix& other) {
            m_buffer = other.m_buffer;
            return *this;
        }
        Matrix& operator=(const std::vector<T>& v) {
            REQUIRE_EQ(v.size(), m_buffer.size(), "cannot assign due to incorrect buffer size");
            m_buffer = v;
            return *this;
        }
        Matrix& operator=(const sparse::Matrix<T>& sparse){
            REQUIRE_GE(m_nrow, sparse.nrow(), "not enough rows in dense matrix ");
            REQUIRE_GT(m_ncol, sparse.max_column_index(), "not enough columns in dense matrix ");
            for (size_t irow = 0ul; irow < sparse.nrow(); ++irow) {
                auto pair = sparse[irow];
                auto nentry = pair.first.size();
                for (size_t ientry = 0ul; ientry < nentry; ++ientry)
                    (*this)(irow, pair.first[ientry]) = pair.second[ientry];
            }
            return *this;
        }

        const size_t& nrow() const {return m_nrow;}
        const size_t& ncol() const {return m_ncol;}

        T* ptr(const size_t &irow=0) {
            return m_buffer.data()+index(irow, 0);
        }

        const T* ptr(const size_t &irow=0) const {
            return m_buffer.data()+index(irow, 0);
        }

        T &operator[](const size_t &iflat) {
            DEBUG_ASSERT_LT(iflat, m_buffer.size(), "index OOB");
            return m_buffer[iflat];
        }

        const T &operator[](const size_t &iflat) const {
            DEBUG_ASSERT_LT(iflat, m_buffer.size(), "index OOB");
            return m_buffer[iflat];
        }

        T &operator()(const size_t &irow, const size_t &icol) {
            return m_buffer[index(irow, icol)];
        }

        const T &operator()(const size_t &irow, const size_t &icol) const {
            return m_buffer[index(irow, icol)];
        }

        void zero() {
            m_buffer.assign(m_buffer.size(), 0);
        }

        bool nearly_equal(const Matrix<T>& other, consts::comp_t<T> eps=consts::eps<T>()) const {
            for (size_t i=0ul; i<m_buffer.size(); ++i) {
                if (!consts::nearly_equal(m_buffer[i], other.m_buffer[i], eps)) return false;
            }
            return true;
        }

    protected:
        /*
         * for bounds safety, we keep these set-from-pointer methods protected
         */
        void set_row(const size_t &irow, const T* v) {
            memcpy(m_buffer.data() + index(irow, 0ul), v, m_ncol * sizeof(T));
        }

        void set_col(const size_t &icol, const T* v) {
            for (size_t irow = 0ul; irow<m_nrow; ++irow) (*this)(irow, icol) = v[irow];
        }

    public:

        /**
         * out-of-place transposition
         */
        void transpose() {
            auto tmp = m_buffer;
            std::swap(m_nrow, m_ncol);
            const auto ptr = tmp.data();
            for (size_t icol = 0ul; icol<m_ncol; ++icol) {
                set_col(icol, ptr);
                ptr+=m_nrow;
            }
        }

        void conj() {
            for (auto& v: m_buffer) v = consts::conj(v);
        }

        void set_row(const size_t &irow, const std::vector<T> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            set_row(irow, v.data());
        }

        void set_col(const size_t &icol, const std::vector<T> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            set_col(icol, v.data());
        }

        EigenSolver<T> diagonalize() const {
            return EigenSolver<T>(*this);
        }

        void mv_multiply(const std::vector<T> &in, std::vector<T> &out);

        std::string to_string() const {
            std::string out;
            for (size_t irow = 0ul; irow < m_nrow; ++irow) {
                for (size_t icol = 0ul; icol < m_ncol; ++icol) {
                    out+=utils::to_string((*this)(irow, icol))+"  ";
                }
                out+="\n";
            }
            return out;
        }
    };

    template<typename T> class SquareMatrix;
    /**
     * double precision, real symmetric eigensolver
     * @param mat
     *  matrix to diagonalize in-place
     * @param evals
     *  vector whose elements hold the eigenvalues
     */
    void diag_inplace(SquareMatrix<double>& mat, std::vector<double>& evals);
    /**
     * double precision, complex hermitian eigensolver
     * @param mat
     *  matrix to diagonalize in-place
     * @param evals
     *  vector whose elements hold the eigenvalues
     */
    void diag_inplace(SquareMatrix<std::complex<double>> &mat, std::vector<double>& evals);

    /**
     * double precision, real non-symmetric eigensolver with no eigenvectors
     * @param mat
     *  matrix to diag_inplace
     * @param evals
     *  vector whose elements hold the eigenvalues
     */
    void diagonalize(SquareMatrix<double>& mat, std::vector<std::complex<double>>& evals);

    template<typename T>
    class SquareMatrix : public Matrix<T> {
    public:
        using Matrix<T>::operator=;
        SquareMatrix(size_t n): Matrix<T>(n, n){}
        SquareMatrix(const sparse::Matrix<T>& sparse) :
            SquareMatrix(std::max(sparse.nrow(), sparse.max_column_index()+1)){
            *this = sparse;
        }

        template<typename U>
        std::vector<U> diag_inplace(){
            std::vector<U> evals;
            dense::diag_inplace(*this, evals);
            return evals;
        }
    };

    /**
     * the storage convention used in this namespace is row-contiguous, but LAPACK assumes column-contiguous ordering.
     * so if we have a product of A (x, y) and B (y, z), we must compute the product
     * BT (z, y) AT (y, x) = CT (z, x)
     */
    struct GemmArgs {
        const char m_transa = 'T', m_transb = 'T';
        const int m_nrow_op_a, m_ncol_op_a, m_ncol_op_b;
        const int& m_nrow_c, &m_ncol_c;
        const int m_lda, m_ldb, m_ldc;

        GemmArgs(int nrowa, int ncola, int nrowb, int ncolb);

        GemmArgs(int nrowa, int ncola, int nrowb, int ncolb, int nrowc, int ncolc);

    };

    void multiply(GemmArgs args, const double *a, const double *b, double *c);

    /**
     * compute the matrix-vector product
     * @param a
     *  matrix with dimensions (m, n)
     * @param b
     *  column vector of length n
     * @param c
     *  column vector of length m
     */
    void multiply(const Matrix<double>& a, const std::vector<double>& b, std::vector<double> &c, char trans='N');

}

#endif //M7_DENSE_H
