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


template<typename T>
class EigenSolver;

extern "C" double ddot_(int *n, double *sx, int *incx, double *sy, int *incy);
extern "C" void
dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
extern "C" void zheev_(char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda, double *w,
                       std::complex<double> *work, int *lwork, double *rwork, int *info);


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

    public:
        const size_t m_nrow, m_ncol;

        Matrix(size_t nrow, size_t ncol) : m_buffer(nrow * ncol, T(0)), m_nrow(nrow), m_ncol(ncol) {}

        Matrix(size_t n) : Matrix(n, n) {}

        T &operator()(const size_t &irow, const size_t &icol) {
            return m_buffer[index(irow, icol)];
        }

        const T &operator()(const size_t &irow, const size_t &icol) const {
            return m_buffer[index(irow, icol)];
        }

        void zero() {
            m_buffer.assign(m_buffer.size(), 0);
        }

        bool is_square() const {
            return m_ncol == m_nrow;
        }

        void set_row(const size_t &irow, const std::vector<T> &v) {
            DEBUG_ASSERT_EQ(v.size(), m_ncol, "length of vector does not match that of matrix row");
            memcpy(m_buffer.data() + index(irow, 0ul), v.data(), m_ncol * sizeof(T));
        }

        EigenSolver<T> diagonalize() const {
            return EigenSolver<T>(*this);
        }

        void multiply(const std::vector<T> &in, std::vector<T> &out);

        void print() {
            std::cout << std::setprecision(10);
            for (size_t irow = 0ul; irow < m_nrow; ++irow) {
                for (size_t icol = 0ul; icol < m_ncol; ++icol) {
                    std::cout << operator()(irow, icol) << "  ";
                }
                std::cout << std::endl;
            }
        }
    };
}

#endif //M7_DENSE_H
