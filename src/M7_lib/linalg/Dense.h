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
#include <M7_lib/util/FpTol.h>
#include <M7_lib/defs.h>
#include <M7_lib/hdf5/Node.h>
#include "M7_lib/linalg/sparse/Dynamic.h"


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

    template<typename T>
    static bool nearly_equal(const T* v1, const T* v2, uint_t size, arith::comp_t<T> rtol, arith::comp_t<T> atol) {
        for (uint_t i = 0ul; i < size; ++i) {
            if (!fptol::near_equal(v1[i], v2[i], rtol, atol)) return false;
        }
        return true;
    }

    template<typename T>
    static bool nearly_equal(const T* v1, const T* v2, uint_t size) {
        return nearly_equal(v1, v2, size, fptol::default_rtol(*v1), fptol::default_atol(*v1));
    }

    template<typename T>
    static bool nearly_equal(const v_t<T> v1, const v_t<T> v2, arith::comp_t<T> rtol, arith::comp_t<T> atol) {
        REQUIRE_EQ(v1.size(), v2.size(), "vectors must have same number of elements to be compared");
        return nearly_equal(v1.data(), v2.data(), v1.size(), rtol, atol);
    }

    template<typename T>
    static bool nearly_equal(const v_t<T> v1, const v_t<T> v2) {
        REQUIRE_EQ(v1.size(), v2.size(), "vectors must have same number of elements to be compared");
        return nearly_equal(v1.data(), v2.data(), v1.size());
    }


    template<typename T>
    class Matrix {
        v_t<T> m_buffer;

        uint_t index(const uint_t &irow, const uint_t &icol) const {
            DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
            DEBUG_ASSERT_LT(icol, m_ncol, "column index OOB");
            return irow * m_ncol + icol;
        }
        /**
         * these can't be const because of the need to swap them in the inplace transpose, so instead they are made
         * private and read-only access is provided by a getter method for each
         */
        uint_t m_nrow, m_ncol;
    public:
        Matrix(uint_t nrow, uint_t ncol) : m_buffer(nrow * ncol, T(0)), m_nrow(nrow), m_ncol(ncol) {
            REQUIRE_TRUE(m_nrow, "matrix must have a non-zero number of rows");
            REQUIRE_TRUE(m_ncol, "matrix must have a non-zero number of columns");
        }

        Matrix(uint_t nrow, const v_t<T>& rows): Matrix(nrow, rows.size()/nrow) {
            *this = rows;
        }

        Matrix(const hdf5::NodeReader& nr, const str_t name) :
            Matrix(nr.dataset_shape(name)[0], nr.dataset_shape(name)[1]) {
            nr.read_data(name, m_buffer.data(), m_buffer.size());
        }
        explicit Matrix(const sparse::dynamic::Matrix<T>& sparse) : Matrix(sparse.nrow(), sparse.max_col_ind() + 1){
            *this = sparse;
        }

        Matrix(const Matrix& other): Matrix(other.m_nrow, other.m_ncol){
            m_buffer = other.m_buffer;
        }
        Matrix& operator=(const Matrix& other) {
            m_buffer = other.m_buffer;
            return *this;
        }

        template<typename U>
        Matrix& operator=(const v_t<U>& v) {
            REQUIRE_EQ(v.size(), m_buffer.size(), "cannot assign due to incorrect buffer size");
            set(v.data());
            return *this;
        }

        Matrix& operator=(const sparse::dynamic::Matrix<T>& sparse){
            REQUIRE_GE(m_nrow, sparse.nrow(), "not enough rows in dense matrix to store contents of source");
            REQUIRE_GT(m_ncol, sparse.max_col_ind(), "not enough columns in dense matrix store contents of source");
            for (uint_t irow = 0ul; irow < sparse.nrow(); ++irow) {
                for (auto& elem : sparse[irow]) (*this)(irow, elem.m_i) = elem.m_v;
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

        const uint_t& nrow() const {return m_nrow;}
        const uint_t& ncol() const {return m_ncol;}
        std::pair<uint_t, uint_t> dims() const {return {nrow(), ncol()};}

        T* ptr(const uint_t &irow=0) {
            return m_buffer.data()+index(irow, 0);
        }

        const T* ptr(const uint_t &irow=0) const {
            return m_buffer.data()+index(irow, 0);
        }

        T &operator[](const uint_t &iflat) {
            DEBUG_ASSERT_LT(iflat, m_buffer.size(), "index OOB");
            return m_buffer[iflat];
        }

        const T &operator[](const uint_t &iflat) const {
            DEBUG_ASSERT_LT(iflat, m_buffer.size(), "index OOB");
            return m_buffer[iflat];
        }

        T &operator()(const uint_t &irow, const uint_t &icol) {
            return m_buffer[index(irow, icol)];
        }

        const T &operator()(const uint_t &irow, const uint_t &icol) const {
            return m_buffer[index(irow, icol)];
        }

        void zero() {
            m_buffer.assign(m_buffer.size(), 0);
        }

        bool nearly_equal(const Matrix<T>& other, arith::comp_t<T> rtol, arith::comp_t<T> atol) const {
            REQUIRE_TRUE(m_buffer.size()==other.m_buffer.size(),
                         "matrices must have same number of elements to be compared");
            return dense::nearly_equal(ptr(), other.ptr(), m_buffer.size(), rtol, atol);
        }

        bool nearly_equal(const Matrix<T>& other) const {
            REQUIRE_TRUE(m_buffer.size()==other.m_buffer.size(),
                         "matrices must have same number of elements to be compared");
            return dense::nearly_equal(ptr(), other.ptr(), m_buffer.size());
        }

        Matrix<T>& operator +=(Matrix<T>& other) {
            DEBUG_ASSERT_EQ(m_nrow, other.m_nrow, "incompatible number of rows");
            DEBUG_ASSERT_EQ(m_ncol, other.m_ncol, "incompatible number of columns");
            for (uint_t ielem = 0ul; ielem < m_buffer.size(); ++ielem) m_buffer[ielem] += other.m_buffer[ielem];
            return *this;
        }

        Matrix<T>& operator -=(Matrix<T>& other) {
            DEBUG_ASSERT_EQ(m_nrow, other.m_nrow, "incompatible number of rows");
            DEBUG_ASSERT_EQ(m_ncol, other.m_ncol, "incompatible number of columns");
            for (uint_t ielem = 0ul; ielem < m_buffer.size(); ++ielem) m_buffer[ielem] -= other.m_buffer[ielem];
            return *this;
        }

    protected:
        /*
         * for bounds safety, we keep these set-from and get-to pointer methods protected
         */


        /**
         * non-converting copy-in
         * @tparam U
         *  U (==T)
         * @param src
         *  data source
         */
        template<typename U>
        typename std::enable_if<std::is_same<U, T>::value, void>::type
        set(const U* src) {
            memcpy(m_buffer.data(), src, m_buffer.size() * sizeof(T));
        }

        /**
         * converting copy-in
         * @tparam U
         *  U (!=T)
         * @param src
         *  data source
         */
        template<typename U>
        typename std::enable_if<!std::is_same<U, T>::value, void>::type
        set_row(const U* src) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            for (uint_t ielem = 0ul; ielem < m_buffer.size(); ++ielem) m_buffer[ielem] = src[ielem];
        }


        /**
         * non-converting copy-in
         * @tparam U
         *  U (==T)
         * @param irow
         *  row index
         * @param src
         *  data source
         */
        template<typename U>
        typename std::enable_if<std::is_same<U, T>::value, void>::type
        set_row(const uint_t &irow, const U* src) {
            memcpy(m_buffer.data() + index(irow, 0ul), src, m_ncol * sizeof(T));
        }

        /**
         * converting copy-in
         * @tparam U
         *  U (!=T)
         * @param irow
         *  row index
         * @param src
         *  data source
         */
        template<typename U>
        typename std::enable_if<!std::is_same<U, T>::value, void>::type
        set_row(const uint_t &irow, const U* src) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            for (uint_t icol = 0ul; icol<m_ncol; ++icol) (*this)(irow, icol) = src[icol];
        }

        /**
         * non-converting copy-out
         * @tparam U
         *  U (==T)
         * @param irow
         *  row index
         * @param dst
         *  data destination
         */
        template<typename U>
        typename std::enable_if<std::is_same<U, T>::value, void>::type
        get_row(const uint_t &irow, U* dst) const {
            memcpy(dst, m_buffer.data() + index(irow, 0ul), m_ncol * sizeof(T));
        }

        /**
         * non-converting copy-out
         * @tparam U
         *  U (!=T)
         * @param irow
         *  row index
         * @param dst
         *  data destination
         */
        template<typename U>
        typename std::enable_if<!std::is_same<U, T>::value, void>::type
        get_row(const uint_t &irow, U* dst) const {
            static_assert(std::is_convertible<T, U>::value, "invalid conversion");
            for (uint_t icol = 0ul; icol<m_ncol; ++icol) dst[icol] = (*this)(irow, icol);
        }

        template<typename U>
        void set_col(const uint_t &icol, const U* src) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            for (uint_t irow = 0ul; irow<m_nrow; ++irow) (*this)(irow, icol) = src[irow];
        }

        template<typename U>
        void get_col(const uint_t &icol, U* dst) const {
            static_assert(std::is_convertible<T, U>::value, "invalid conversion");
            for (uint_t irow = 0ul; irow<m_nrow; ++irow) dst[irow] = (*this)(irow, icol);
        }

    public:

        /**
         * in-place physical transposition via out-of-place buffer procedure
         */
        void transpose() {
            auto tmp = m_buffer;
            std::swap(m_nrow, m_ncol);
            auto ptr = tmp.data();
            for (uint_t icol = 0ul; icol<m_ncol; ++icol) {
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

        template<typename U>
        void set_row(const uint_t &irow, const v_t<U> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            // copies byte wise if U==T, else converts element wise
            set_row(irow, v.data());
        }

        template<typename U>
        void set_col(const uint_t &icol, const v_t<U> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            // copies element wise if U==T, else converts element wise
            set_col(icol, v.data());
        }

        template<typename U>
        void get_row(const uint_t &irow, v_t<U> &v) const {
            v.resize(m_ncol);
            // copies byte wise if U==T, else converts element wise
            get_row(irow, v.data());
        }

        template<typename U>
        void get_col(const uint_t &icol, v_t<U> &v) const {
            v.resize(m_nrow);
            // copies element wise if U==T, else converts element wise
            get_col(icol, v.data());
        }

        str_t to_string() const {
            str_t out;
            for (uint_t irow = 0ul; irow < m_nrow; ++irow) {
                for (uint_t icol = 0ul; icol < m_ncol; ++icol) {
                    out+=convert::to_string((*this)(irow, icol))+"  ";
                }
                out+="\n";
            }
            return out;
        }

        void save(str_t name, hdf5::NodeWriter& nw, uint_t irank=0) const {
            uintv_t shape;
            shape.push_back(m_nrow);
            shape.push_back(m_ncol);
            nw.write_data(name, m_buffer, shape, {"nrecord", "ncol"}, irank);
        }
    };

    template<typename T>
    class SquareMatrix : public Matrix<T> {
        static uint_t nrow_from_flat_size(uint_t n) {
            const auto nrow = integer::sqrt(n);
            REQUIRE_NE(nrow, ~0ul, "SquareMatrix initialization from std::vector requires square number of elements");
            return nrow;
        }
    public:
        using Matrix<T>::operator=;
        SquareMatrix(uint_t n): Matrix<T>(n, n){}
        SquareMatrix(const v_t<T>& v): Matrix<T>(nrow_from_flat_size(v.size()), v){}
        SquareMatrix(const sparse::dynamic::Matrix<T>& sparse) :
                SquareMatrix(std::max(sparse.nrow(), sparse.max_col_ind()+1)){
            *this = sparse;
        }
        SquareMatrix(const SquareMatrix& other) : SquareMatrix(other.nrow()){
            *this = other;
        }
        SquareMatrix& operator=(const SquareMatrix& other) {
            Matrix<T>::operator=(other);
            return *this;
        }

        void symmetrize(bool conj=dtype::is_complex<T>()) {
            const auto n = Matrix<T>::nrow();
            for (uint_t irow=0ul; irow < n; ++irow) {
                for (uint_t icol=0ul; icol < n; ++icol) {
                    auto this_element = (*this)(irow, icol);
                    if (!this_element) continue;
                    if (conj) this_element = arith::conj(this_element);
                    auto& other_element = (*this)(icol, irow);
                    if (other_element) {
                        DEBUG_ASSERT_NEARLY_EQ(this_element, other_element, "inconsistent elements");
                    }
                    if (irow!=icol) other_element = this_element;
                }
            }
        }

        bool is_diagonal() const {
            const auto n = Matrix<T>::ncol();
            if (n==1) return true;
            auto first_row_ptr = Matrix<T>::ptr()+1ul;
            if (std::any_of(first_row_ptr, first_row_ptr+n, [](const T& v){return v;})) return false;
            /*
             * non-diagonal part of first row and first element of second row are all zero, so compare all subsequent
             * non-diagonal parts with this string
             */
            auto ptr = first_row_ptr;
            for (auto irow=2ul; irow < n; ++irow) {
                ptr+=n+1;
                if (std::memcmp(ptr, first_row_ptr, n*sizeof(T))) return false;
            }
            return true;
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
        // nrecord of op(A)
        const int m_m;
        // ncol of op(B)
        const int m_n;
        // ncol of op(A)
        const int m_k;
        // nrecord of op(B)
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
        GemmWrapper(uint_t nrowp, uint_t ncolp, uint_t nrowq, uint_t ncolq, char transp, char transq, uint_t nrowr, uint_t ncolr);

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
    void multiply(const Matrix<T>& p, const v_t<T>& q, v_t<T>& r,
                  char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        r.resize(transp=='N' ? p.nrow() : p.ncol());
        GemmWrapper wrapper(p.nrow(), p.ncol(), q.size(), 1ul, transp, transq, r.size(), 1ul);
        wrapper.multiply(p.ptr(), q.data(), r.data(), alpha, beta);
    }

    template<typename T>
    v_t<T> multiply(const Matrix<T>& p, const v_t<T>& q,
                                char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        v_t<T> r;
        multiply(p, q, r, transp, transq, alpha, beta);
        return r;
    }

    template<typename T>
    T inner_product(const v_t<T>& p, const v_t<T>& q, bool herm=false){
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

    bool diag(const SquareMatrix<float>& mat, v_t<float>& evals);

    bool diag(const SquareMatrix<float>& mat, SquareMatrix<float>& evecs, v_t<float>& evals);

    bool diag(const SquareMatrix<float>& mat, v_t<std::complex<float>>& evals);


    bool diag(const SquareMatrix<double>& mat, v_t<double>& evals);

    bool diag(const SquareMatrix<double>& mat, SquareMatrix<double>& evecs, v_t<double>& evals);

    bool diag(const SquareMatrix<double>& mat, v_t<std::complex<double>>& evals);



    bool diag(const SquareMatrix<std::complex<double>>& mat, v_t<double>& evals);


    bool diag(const SquareMatrix<std::complex<double>>& mat,
              SquareMatrix<std::complex<double>>& evecs, v_t<double>& evals);

    bool diag(const SquareMatrix<std::complex<double>>& mat, v_t<double>& evals);
}

#endif //M7_DENSE_H
