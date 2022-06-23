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

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/consts.h>
#include <M7_lib/defs.h>
#include <M7_lib/io/HDF5Wrapper.h>
#include "Sparse.h"


extern "C" void ssyev_(const char *jobz, const char *uplo, const int *n, float *a, const int *lda, float *w,
                       float *work, const int *lwork, int *info);

extern "C" void sgeev_(const char *jobvl, const char *jobvr, const int *n, float *a, const int *lda,
                       float *wr, float *wi, float *vl, const int *ldvl, float *vr, const int *ldvr,
                       float *work, const int *lwork, int *info);


extern "C" void dsyev_(const char *jobz, const char *uplo, const int *n, double *a, const int *lda, double *w,
                       double *work, const int *lwork, int *info);

extern "C" void zheev_(const char *jobz, const char *uplo, const int *n, std::complex<double> *a, const int *lda,
                       double *w, std::complex<double> *work, const int *lwork, double *rwork, int *info);


extern "C" void dgeev_(const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda,
                       double *wr, double *wi,
                       double *vl, const int* ldvl,
                       double *vr, const int* ldvr,
                       double *work, const int *lwork, int *info);





extern "C" void sgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const float *alpha,
                       const float *a, const int *lda,
                       const float *b, const int *ldb,
                       const float *beta,
                       float *c, const int *ldc);

extern "C" void dgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const double *alpha,
                       const double *a, const int *lda,
                       const double *b, const int *ldb,
                       const double *beta,
                       double *c, const int *ldc);

extern "C" void cgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const std::complex<float> *alpha,
                       const std::complex<float> *a, const int *lda,
                       const std::complex<float> *b, const int *ldb,
                       const std::complex<float> *beta,
                       std::complex<float> *c, const int *ldc);

extern "C" void zgemm_(const char *transa, const char *transb,
                       const int *m, const int *n, const int *k,
                       const std::complex<double> *alpha,
                       const std::complex<double> *a, const int *lda,
                       const std::complex<double> *b, const int *ldb,
                       const std::complex<double> *beta,
                       std::complex<double> *c, const int *ldc);


namespace dense {

    using namespace fptol;
    template<typename T>
    static bool nearly_equal(const T* v1, const T* v2, size_t size, comp_t<T> atol = default_atol_near<T>()) {
        for (size_t i = 0ul; i < size; ++i) {
            if (!fptol::nearly_equal(v1[i], v2[i], comp_t<T>(0), atol)) return false;
        }
        return true;
    }

    template<typename T>
    static bool nearly_equal(const std::vector<T> v1, const std::vector<T> v2, comp_t<T> atol = default_atol_near<T>()) {
        REQUIRE_EQ(v1.size(), v2.size(), "vectors must have same number of elements to be compared");
        return nearly_equal(v1.data(), v2.data(), v1.size(), atol);
    }


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

        /**
         * byte-wise equality
         * @param other
         *  matrix to compare against
         * @return
         *  true only if both buffers are identical
         */
        bool operator==(const Matrix<T>& other) const {
            return !std::memcmp(ptr(), other.ptr(), sizeof(T)*m_ncol*m_nrow);
        }

        const size_t& nrow() const {return m_nrow;}
        const size_t& ncol() const {return m_ncol;}
        std::pair<size_t, size_t> dims() const {return {nrow(), ncol()};}

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

        bool nearly_equal(const Matrix<T>& other, arith::comp_t<T> atol= fptol::default_atol_near<T>()) const {
            REQUIRE_TRUE(m_buffer.size()==other.m_buffer.size(),
                         "matrices must have same number of elements to be compared");
            return dense::nearly_equal(ptr(), other.ptr(), m_buffer.size(), atol);
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
         * in-place physical transposition via out-of-place buffer procedure
         */
        void transpose() {
            auto tmp = m_buffer;
            std::swap(m_nrow, m_ncol);
            auto ptr = tmp.data();
            for (size_t icol = 0ul; icol<m_ncol; ++icol) {
                set_col(icol, ptr);
                ptr+=m_nrow;
            }
        }

        /**
         * in-place complex conjugation
         */
        void conj() {
            for (auto& v: m_buffer) v = arith::conj(v);
        }

        /**
         * in-place hermitian conjugate-transposition
         */
        void dagger() {
            transpose();
            conj();
        }

        void set_row(const size_t &irow, const std::vector<T> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            set_row(irow, v.data());
        }

        void set_col(const size_t &icol, const std::vector<T> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            set_col(icol, v.data());
        }

        std::string to_string() const {
            std::string out;
            for (size_t irow = 0ul; irow < m_nrow; ++irow) {
                for (size_t icol = 0ul; icol < m_ncol; ++icol) {
                    out+=utils::convert::to_string((*this)(irow, icol))+"  ";
                }
                out+="\n";
            }
            return out;
        }

        void save(std::string name, hdf5::GroupWriter& gw, size_t irank=0) const {
            defs::inds_t shape;
            shape.push_back(m_nrow);
            shape.push_back(m_ncol);
            gw.save(name, m_buffer, shape, {"nrow", "ncol"}, irank);
        }
    };

    template<typename T>
    class SquareMatrix : public Matrix<T> {
    public:
        using Matrix<T>::operator=;
        SquareMatrix(size_t n): Matrix<T>(n, n){}
        SquareMatrix(const sparse::Matrix<T>& sparse) :
                SquareMatrix(std::max(sparse.nrow(), sparse.max_column_index()+1)){
            *this = sparse;
        }
        SquareMatrix(const SquareMatrix& other) : SquareMatrix(other.nrow()){
            *this = other;
        }
        SquareMatrix& operator=(const SquareMatrix& other) {
            Matrix<T>::operator=(other);
            return *this;
        }
    };


    /**
     * Each of the overloaded GEMM implementations has a common interface in terms of dimensional extents, strides,
     * and transposition. These are expressed in the members of this class to reduce code duplication, and to handle
     * the required reversal of the order of multiplication.
     *
     * BLAS assumes the Fortran convention of column-major ordering, i.e. columns are assumed to be consecutively laid-
     * out in linear memory. This is in contrast to the C/C++ convention of row-major ordering. Where we have a (m x n)
     * matrix in row-major memory, we have a column-major representation of its transpose. Thus, row and column labels
     * are swapped as well as the identities of the buffers containing the matrices before passing to GEMM.
     *
     * let the desired product of matrices and the GEMM product be denoted by P.Q = R and A.B = C respectively.
     * we actually do QT.PT = RT with GEMM.
     */
    struct GemmWrapper {
        const char m_transa;
        const char m_transb;
        // nrow of op(A)
        const int m_m;
        // ncol of op(B)
        const int m_n;
        // ncol of op(A)
        const int m_k;
        // nrow of op(B)
        const int m_nrow_opb;
        // strides of row-contiguous arrays
        const int m_lda, m_ldb, m_ldc;

        /**
         * @param t
         *  input transposition flag for validation
         * @return
         *  the capitalized transposition flag
         */
        static char valid_trans(char t);

        /**
         * overloaded ctor to conveniently check dimensional compatibility of result
         */
        GemmWrapper(size_t nrowp, size_t ncolp, size_t nrowq, size_t ncolq, char transp, char transq, size_t nrowr, size_t ncolr);

        void multiply(const float* p, const float* q, float* r, float alpha, float beta);

        void multiply(const double* p, const double* q, double* r, double alpha, double beta);

        void multiply(const std::complex<float>* p, const std::complex<float>* q, std::complex<float>* r,
                      std::complex<float> alpha, std::complex<float> beta);

        void multiply(const std::complex<double>* p, const std::complex<double>* q, std::complex<double>* r,
                      std::complex<double> alpha, std::complex<double> beta);


    };

    template<typename T>
    void multiply(const Matrix<T>& p, const Matrix<T>& q, Matrix<T>& r,
                  char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        GemmWrapper wrapper(p.nrow(), p.ncol(), q.nrow(), q.ncol(), transp, transq, r.nrow(), r.ncol());
        wrapper.multiply(p.ptr(), q.ptr(), r.ptr(), alpha, beta);
    }

    template<typename T>
    void multiply(const Matrix<T>& p, const std::vector<T>& q, std::vector<T>& r,
                  char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        r.resize(transp=='N' ? p.nrow() : p.ncol());
        GemmWrapper wrapper(p.nrow(), p.ncol(), q.size(), 1ul, transp, transq, r.size(), 1ul);
        wrapper.multiply(p.ptr(), q.data(), r.data(), alpha, beta);
    }

    template<typename T>
    std::vector<T> multiply(const Matrix<T>& p, const std::vector<T>& q,
                                char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        std::vector<T> r;
        multiply(p, q, r, transp, transq, alpha, beta);
        return r;
    }

    template<typename T>
    T inner_product(const std::vector<T>& p, const std::vector<T>& q, bool herm=false){
        REQUIRE_EQ(p.size(), q.size(), "vector lengths don't match");
        GemmWrapper wrapper(p.size(), 1ul, q.size(), 1ul, herm ? 'C' : 'T', 'N', 1ul, 1ul);
        T r = 0;
        wrapper.multiply(p.data(), q.data(), &r, 1.0, 0.0);
        return r;
    }

    /*
     * Diagonalization routines do not admit the same common interface as GEMM, so here there is no wrapper class, but
     * we use overloading to infer the LAPACK routine which should be called in order to perform the diagonalization.
     * All of these LAPACK routines are destructive of the input matrix, but all diagonalization functions defined in
     * this namespace make a copy to pass to the solver and so do not modify the given matrix. Functions are not
     * provided for the complex-symmetric case, only the hermitian case.
     *
     * diag methods return true if LAPACK call was successful
     *
     *    A (matrix)      R (right evecs)  L (left evecs)  D (evals)       Routine
     *    float           -                -               float           ssyev
     *    float           float            -               float           ssyev
     *    float           -                -               complex float   sgeev
     *    float           complex float    complex float   complex float   sgeev
     *
     *    double          -                -               double          dsyev
     *    double          double           -               double          dsyev
     *    double          -                -               complex double  dgeev
     *    double          complex double   complex double  complex double  dgeev
     *
     *    complex float   -                -               float           cheev
     *    complex float   complex float    -               float           cheev
     *    complex float   -                -               complex float   cgeev
     *    complex float   complex float    complex float   complex float   cgeev
     *
     *    complex double  -                -               double          zheev
     *    complex double  complex double   -               double          zheev
     *    complex double  -                -               complex double  zgeev
     *    complex double  complex double   complex double  complex double  zgeev
     */

    bool diag(const SquareMatrix<float>& mat, std::vector<float>& evals);

    bool diag(const SquareMatrix<float>& mat, SquareMatrix<float>& evecs, std::vector<float>& evals);

    bool diag(const SquareMatrix<float>& mat, std::vector<std::complex<float>>& evals);


    bool diag(const SquareMatrix<double>& mat, std::vector<double>& evals);

    bool diag(const SquareMatrix<double>& mat, SquareMatrix<double>& evecs, std::vector<double>& evals);

    bool diag(const SquareMatrix<double>& mat, std::vector<std::complex<double>>& evals);



    bool diag(const SquareMatrix<std::complex<double>>& mat, std::vector<double>& evals);


    bool diag(const SquareMatrix<std::complex<double>>& mat,
              SquareMatrix<std::complex<double>>& evecs, std::vector<double>& evals);

    bool diag(const SquareMatrix<std::complex<double>>& mat, std::vector<double>& evals);
}

#endif //M7_DENSE_H
