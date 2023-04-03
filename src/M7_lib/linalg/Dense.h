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
#include <M7_lib/table/Buffer.h>
#include <M7_lib/util/Sort.h>
#include "M7_lib/linalg/sparse/Dynamic.h"
#include "M7_lib/hdf5/Node.h"
#include "M7_lib/hdf5/DatasetTransaction.h"


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

extern "C" void zgeev_(const char* jobvl, const char* jobvr, const int* n, std::complex<double>* a,
                       const int* lda, std::complex<double>* w, std::complex<double>* vl, const int* ldvl,
                       std::complex<double>* vr, const int* ldvr, std::complex<double>* work,
                       const int* lwork, double* rwork, int* info);



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


    class MatrixBase {
        Buffer m_buffer;
    protected:
        Buffer::Window m_bw;
        const uint_t m_element_size;
        /**
         * these can't be const because of the need to swap them in the inplace transpose, so instead they are made
         * private and read-only access is provided by a getter method for each
         */
        uint_t m_nrow = 0ul;
        uint_t m_ncol = 0ul;
        uint_t m_nelement = 0ul;
        uint_t m_row_size = 0ul;
        // size is not necessarily m_bw.m_size, since buffer windows are padded to a whole number of words
        uint_t m_size = 0ul;

        /**
         * obviously, in non-node_sharing mode, an MPI rank has the authority to modify all elements
         * but otherwise, we still allow element-wise modification by non-root-on-node ranks, but global modification
         * e.g. set to zero, scale, add/subtract, should all take place on the node root
         * @return
         *  true if this rank participates in globally-modifying operations of the buffer
         */
        bool i_can_globally_modify() const {
            return !m_buffer.m_node_shared || mpi::on_node_i_am_root();
        }

        void set_sizes(uint_t nrow, uint_t ncol);

        buf_t* begin() {
            return m_bw.begin();
        }
        const buf_t* cbegin() const {
            return m_bw.cbegin();
        }

        buf_t* begin(uint_t irow) {
            DEBUG_ASSERT_LT(irow, m_nrow, "index OOB");
            return begin() + irow * m_row_size;
        }
        const buf_t* cbegin(uint_t irow) const {
            DEBUG_ASSERT_LT(irow, m_nrow, "index OOB");
            return cbegin() + irow * m_row_size;
        }

        buf_t* begin(uint_t irow, uint_t icol) {
            DEBUG_ASSERT_LT(icol, m_ncol, "index OOB");
            return begin(irow) + icol * m_element_size;
        }
        const buf_t* cbegin(uint_t irow, uint_t icol) const {
            DEBUG_ASSERT_LT(icol, m_ncol, "index OOB");
            return cbegin(irow) + icol * m_element_size;
        }

        template<typename T>
        T* begin_as() {
            return reinterpret_cast<T*>(begin());
        }

        template<typename T>
        const T* cbegin_as() const {
            return reinterpret_cast<const T*>(cbegin());
        }

        template<typename T>
        T* begin_as(uint_t irow) {
            return reinterpret_cast<T*>(begin(irow));
        }

        template<typename T>
        const T* cbegin_as(uint_t irow) const {
            return reinterpret_cast<const T*>(cbegin(irow));
        }

    public:

        bool compatible(const MatrixBase& other) const;

        void resize(uint_t nrow, uint_t ncol);

        MatrixBase(uint_t nrow, uint_t ncol, uint_t element_size, bool node_shared);

        MatrixBase(const MatrixBase& other);

        MatrixBase& operator=(const MatrixBase& other);

        /**
         * byte-wise equality
         * @param other
         *  matrix to compare against
         * @return
         *  true only if both buffers are identical
         */
        bool operator==(const MatrixBase& other) const;

        uint_t nrow() const {return m_nrow;}
        uint_t ncol() const {return m_ncol;}

        void zero();

        /**
         * in-place physical transposition via out-of-place buffer procedure
         */
        void transpose();

        void reorder_rows(const uintv_t& order);

        /**
         * byte-wise getters and setters are kept protected. subclass methods with std::vector<T> args are safer public
         * interfaces
         */
    protected:
        void set(const void* src);

        void set_row(uint_t irow, const void* src);

        void get_row(uint_t irow, void* dst) const;

        void set_col(uint_t icol, const void* src);

        void get_col(uint_t icol, void* dst);
    };

    template<typename T>
    class Matrix : public MatrixBase {

    protected:

        T* tbegin(uint_t irow) {
            return begin_as<T>(irow);
        }

        const T* ctbegin(uint_t irow) const {
            return cbegin_as<T>(irow);
        }

    public:
        T* tbegin() {
            return begin_as<T>();
        }

        const T* ctbegin() const {
            return cbegin_as<T>();
        }

        Matrix(uint_t nrow, uint_t ncol, bool node_shared=false):
            MatrixBase(nrow, ncol, sizeof(T), node_shared){}

        Matrix(uint_t nrow, const v_t<T>& rows): Matrix(nrow, rows.size()/nrow) {
            *this = rows;
        }

        Matrix(const hdf5::NodeReader& nr, const str_t name, bool this_rank, bool node_shared=false):
                Matrix(0, 0, node_shared){
            hdf5::DatasetLoader dl(nr, name, false, this_rank);
            const auto& shape = dl.m_format.m_h5_shape;
            resize(shape[0], shape[1]);
            std::list<hdf5::Attr> attrs;
            dl.load_dist_list(nr, name, begin(), m_size, hdf5::Type::make<T>(), dtype::is_complex<T>(),
                false, this_rank, attrs);
        }

        explicit Matrix(const sparse::dynamic::Matrix<T>& sparse) : Matrix(sparse.nrow(), sparse.max_col_ind() + 1){
            *this = sparse;
        }

        template<typename U>
        Matrix& operator=(const v_t<U>& v) {
            REQUIRE_EQ(v.size(), m_nelement, "cannot assign due to incorrect buffer size");
            set(v.data());
            return *this;
        }

        Matrix& operator=(const sparse::dynamic::Matrix<T>& sparse){
            const auto nrow = std::max(m_nrow, sparse.nrow());
            const auto ncol = std::max(m_ncol, sparse.max_col_ind()+1);
            if ((nrow != m_nrow) || (ncol != m_ncol)) resize(nrow, ncol);
            if (i_can_globally_modify()) {
                for (uint_t irow = 0ul; irow < sparse.nrow(); ++irow) {
                    for (auto &elem: sparse[irow]) (*this)(irow, elem.m_i) = elem.m_v;
                }
            }
            return *this;
        }

        T &operator[](uint_t ielement) {
            DEBUG_ASSERT_LT(ielement, m_nelement, "element index OOB");
            return tbegin()[ielement];
        }

        const T &operator[](uint_t ielement) const {
            DEBUG_ASSERT_LT(ielement, m_nelement, "element index OOB");
            return ctbegin()[ielement];
        }

        T &operator()(uint_t irow, uint_t icol) {
            return *reinterpret_cast<T*>(begin(irow, icol));
        }

        const T &operator()(uint_t irow, uint_t icol) const {
            return *reinterpret_cast<const T*>(cbegin(irow, icol));
        }

        bool nearly_equal(const Matrix<T>& other, arith::comp_t<T> rtol, arith::comp_t<T> atol) const {
            REQUIRE_TRUE(compatible(other), "matrices must be compatible to be compared");
            return dense::nearly_equal(ctbegin(), other.ctbegin(), m_nelement, rtol, atol);
        }

        bool nearly_equal(const Matrix<T>& other) const {
            REQUIRE_TRUE(compatible(other), "matrices must be compatible to be compared");
            return dense::nearly_equal(ctbegin(), other.ctbegin(), m_nelement);
        }

        Matrix<T>& operator +=(const Matrix<T>& other) {
            if (i_can_globally_modify()) {
                REQUIRE_TRUE(compatible(other), "matrices must be compatible to be added");
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) (*this)[ielem] += other[other];
            }
            return *this;
        }

        Matrix<T>& operator -=(const Matrix<T>& other) {
            if (i_can_globally_modify()) {
                REQUIRE_TRUE(compatible(other), "matrices must be compatible to be added");
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) (*this)[ielem] -= other[other];
            }
            return *this;
        }

        Matrix<T>& operator *=(T scalar) {
            if (i_can_globally_modify()) {
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) (*this)[ielem] *= scalar;
            }
            return *this;
        }

        Matrix<T>& operator /=(T scalar) {
            if (i_can_globally_modify()) {
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) (*this)[ielem] /= scalar;
            }
            return *this;
        }


    protected:
        /*
         * for bounds safety, we keep these set-from and get-to pointer methods protected
         */

        /**
         * non-converting copy-in
         * @param src
         *  data source
         */
        void set(const T* src, tag::Int<0> /*convert*/) {
            MatrixBase::set(src);
        }

        /**
         * converting copy-in
         * @tparam U
         *  U (!=T)
         * @param src
         *  data source
         */
        template<typename U>
        void set(const U* src, tag::Int<1> /*convert*/) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            if (i_can_globally_modify()) {
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) (*this)[ielem] = src[ielem];
            }
        }

        /*
         * dispatch the appropriate tagged method
         */
        template<typename U>
        void set(const U* src) {
            set(src, tag::Int<!std::is_same<U, T>::value>());
        }

        /**
         * non-converting copy-in
         * @param irow
         *  row index
         * @param src
         *  data source
         */
        void set_row(uint_t irow, const T* src, tag::Int<0> /*convert*/) {
            MatrixBase::set_row(irow, src);
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
        void set_row(uint_t irow, const U* src, tag::Int<1> /*convert*/) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            if (i_can_globally_modify()) {
                for (uint_t icol = 0ul; icol < m_ncol; ++icol) (*this)(irow, icol) = src[icol];
            }
        }

        /*
         * dispatch the appropriate tagged method
         */
        template<typename U>
        void set_row(uint_t irow, const U* src) {
            set_row(irow, src, tag::Int<!std::is_same<U, T>::value>());
        }

        /**
         * non-converting copy-out
         * @param irow
         *  row index
         * @param dst
         *  data destination
         */
        void get_row(uint_t irow, T* dst, tag::Int<0> /*convert*/) const {
            MatrixBase::get_row(irow, dst);
        }

        /**
         * converting copy-out
         * @tparam U
         *  U (!=T)
         * @param irow
         *  row index
         * @param dst
         *  data destination
         */
        template<typename U>
        void get_row(uint_t irow, U* dst, tag::Int<1> /*convert*/) const {
            static_assert(std::is_convertible<T, U>::value, "invalid conversion");
            for (uint_t icol = 0ul; icol < m_ncol; ++icol) dst[icol] = (*this)(irow, icol);
        }

        /*
         * dispatch the appropriate tagged method
         */
        template<typename U>
        void get_row(uint_t irow, U* dst) const {
            get_row(irow, dst, tag::Int<!std::is_same<U, T>::value>());
        }

        /*
         * non-contiguous columns means there is no point in differentiating between converting and non-converting ops
         */
        template<typename U>
        void set_col(uint_t icol, const U* src) {
            static_assert(std::is_convertible<U, T>::value, "invalid conversion");
            if (i_can_globally_modify()) {
                for (uint_t irow = 0ul; irow < m_nrow; ++irow) (*this)(irow, icol) = src[irow];
            }
        }

        template<typename U>
        void get_col(uint_t icol, U* dst) const {
            static_assert(std::is_convertible<T, U>::value, "invalid conversion");
            if (i_can_globally_modify()) {
                for (uint_t irow = 0ul; irow < m_nrow; ++irow) dst[irow] = (*this)(irow, icol);
            }
        }

    public:

        /**
         * in-place complex conjugation
         */
        void conj() {
            if (i_can_globally_modify()) {
                for (auto it = tbegin(); it != tbegin() + m_nelement; ++it) *it = arith::conj(*it);
            }
        }

        /**
         * set every element with magnitude less than thresh to zero
         */
        void screen(const arith::comp_t<T>& thresh) {
            if (i_can_globally_modify()) {
                for (uint_t ielem = 0ul; ielem < m_nelement; ++ielem) {
                    auto& elem = (*this)[ielem];
                    if (std::abs(elem) < thresh) elem = 0.0;
                }
            }
        }

        /**
         * in-place hermitian conjugate-transposition
         */
        void dagger() {
            transpose();
            conj();
        }

        void all_sum() {
            // todo: make node-sharing safe
            Matrix<T> tmp = *this;
            mpi::all_sum(tmp.ctbegin(), tbegin(), m_nelement);
        }

        template<typename U>
        void set_row(uint_t irow, const v_t<U> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            // copies byte wise if U==T, else converts element wise
            set_row(irow, v.data());
        }

        template<typename U>
        void set_col(uint_t icol, const v_t<U> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            // copies element wise if U==T, else converts element wise
            set_col(icol, v.data());
        }

        template<typename U>
        void get_row(uint_t irow, v_t<U> &v) const {
            v.resize(m_ncol);
            get_row(irow, v.data());
        }

        template<typename U>
        void get_col(uint_t icol, v_t<U> &v) const {
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
            hdf5::DatasetSaver::save_array(nw, name, ctbegin(), shape, {"row", "col"}, irank);
        }
    };

    template<typename T>
    class Vector : public Matrix<T> {
        using MatrixBase::m_bw;
        using MatrixBase::m_ncol;
        using MatrixBase::m_nrow;
        using MatrixBase::set_sizes;
        using MatrixBase::i_can_globally_modify;
    public:
        explicit Vector(uint_t nelement, bool node_shared=false): Matrix<T>(1, nelement, node_shared){}

        explicit Vector(const v_t<T>& v): Vector(v.size()) {
            Matrix<T>::operator=(v);
        }

        uint_t nelement() const {
            return MatrixBase::m_nelement;
        }

        void reorder(const uintv_t& order) {
            if (i_can_globally_modify())
                sort::reorder(MatrixBase::begin(), MatrixBase::m_element_size, order);
            if (m_bw.node_shared()) mpi::barrier_on_node();
        }

        void sort_inplace(bool asc, bool absval) {
            if (i_can_globally_modify())
                sort::inplace(Matrix<T>::tbegin(), MatrixBase::m_nelement, asc, absval);
            if (m_bw.node_shared()) mpi::barrier_on_node();
        }

        Vector<T>& sorted(bool asc, bool absval) {
            if (i_can_globally_modify())
                sort::inplace(Matrix<T>::tbegin(), MatrixBase::m_nelement, asc, absval);
            if (m_bw.node_shared()) mpi::barrier_on_node();
            return *this;
        }

        uintv_t sort_inds(bool asc, bool absval) const {
            return sort::inds(Matrix<T>::ctbegin(), MatrixBase::m_nelement, asc, absval);
        }

        void sort_inds(uintv_t& order, bool asc, bool absval) const {
            sort::inds(order, Matrix<T>::ctbegin(), MatrixBase::m_nelement, asc, absval);
        }

        Vector(const hdf5::NodeReader& nr, const str_t name, bool this_rank, bool node_shared=false):
                Matrix<T>(nr, name, this_rank, node_shared){
            // vectors must have one row
            if (m_nrow!=1) set_sizes(m_ncol, m_nrow);
            REQUIRE_EQ(m_nrow, 1ul, "read data has non-vector shape");
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

        void resize(uint_t n) {
            Matrix<T>::resize(n, n);
        }

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

        bool is_symmetric(bool conj=dtype::is_complex<T>()) {
            const auto n = Matrix<T>::nrow();
            for (uint_t irow=0ul; irow < n; ++irow) {
                const auto& diag = (*this)(irow, irow);
                if (conj && arith::imag(diag)!=0.0) return false;
                for (uint_t icol = irow + 1; icol < n; ++icol) {
                    const auto& upper = (*this)(irow, icol);
                    const auto& lower = (*this)(icol, irow);
                    if (conj) {
                        if (upper != arith::conj(lower)) return false;
                    }
                    else {
                        if (upper != lower) return false;
                    }
                }
            }
            // no non-symmetric pairs found
            return true;
        }

        void symmetrize(bool conj=dtype::is_complex<T>()) {
            if (!MatrixBase::i_can_globally_modify()) return;
            const auto n = Matrix<T>::nrow();
            for (uint_t irow=0ul; irow < n; ++irow) {
                for (uint_t icol=0ul; icol < n; ++icol) {
                    auto this_element = (*this)(irow, icol);
                    if (this_element == 0.0) continue;
                    if (conj) this_element = arith::conj(this_element);
                    auto& other_element = (*this)(icol, irow);
                    if (other_element != 0.0) {
                        DEBUG_ASSERT_NEAR_EQ(this_element, other_element, "inconsistent elements");
                    }
                    if (irow!=icol) other_element = this_element;
                }
            }
        }

        bool is_diagonal() const {
            const auto n = Matrix<T>::ncol();
            if (n==1) return true;
            auto first_row_ptr = Matrix<T>::ctbegin()+1ul;
            if (std::any_of(first_row_ptr, first_row_ptr+n, [](const T& v){return v!=0.0;})) return false;
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

        v_t<T> get_diagonal() const {
            const auto n = Matrix<T>::ncol();
            v_t<T> tmp;
            tmp.reserve(n);
            for (uint_t i=0ul; i<n; ++i) tmp.push_back((*this)(i, i));
            return tmp;
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
        wrapper.multiply(p.ctbegin(), q.ctbegin(), r.tbegin(), alpha, beta);
    }

    template<typename T>
    void multiply(const Matrix<T>& p, const v_t<T>& q, v_t<T>& r,
                  char transp='N', char transq='N', T alpha=1.0, T beta=0.0){
        r.resize(transp=='N' ? p.nrow() : p.ncol());
        GemmWrapper wrapper(p.nrow(), p.ncol(), q.size(), 1ul, transp, transq, r.size(), 1ul);
        wrapper.multiply(p.ctbegin(), q.data(), r.data(), alpha, beta);
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


    bool diag(const dense::SquareMatrix<std::complex<double>> &mat, v_t<std::complex<double>> &evals);

}

#endif //M7_DENSE_H