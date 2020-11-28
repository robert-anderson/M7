//
// Created by Robert John Anderson on 2020-01-24.
//

#ifndef M7_MATRIX_H
#define M7_MATRIX_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cstring>
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

template<typename T>
class Matrix {
    std::vector<T> m_data;

    T *data(const size_t &irow, const size_t &icol) {
        ASSERT(icol * m_nrow + irow < m_data.size());
        return m_data.data() + icol * m_nrow + irow;
    }

public:
    const size_t m_nrow, m_ncol;

    Matrix(size_t nrow, size_t ncol) : m_data(nrow * ncol, T(0)), m_nrow(nrow), m_ncol(ncol) {}

    Matrix(size_t n) : Matrix(n, n) {}

    T &operator()(const size_t &irow, const size_t &icol) {
        ASSERT(irow < m_nrow);
        ASSERT(icol < m_ncol);
        return *data(irow, icol);
    }

    void zero(){
        m_data.assign(m_data.size(), 0);
    }

    bool is_square() const {
        return m_ncol == m_nrow;
    }

    void set_row(const size_t &irow, const std::vector<T> &v) {
        ASSERT(v.size() == m_ncol);
        memcpy(data(irow), v.data(), m_ncol * sizeof(T));
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


#endif //M7_MATRIX_H
