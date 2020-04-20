//
// Created by Robert John Anderson on 2020-01-24.
//

#ifndef M7_MATRIX_H
#define M7_MATRIX_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include "src/consts.h"
#include "src/defs.h"

template<typename T>
class EigenSolver;

extern "C" double ddot_(int *n, double *sx, int *incx, double *sy, int *incy);
extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
extern "C" void zheev_( char* jobz, char* uplo, int* n, std::complex<double>* a, int* lda, double* w,
                        std::complex<double>* work, int* lwork, double* rwork, int* info);

template<typename T>
class Matrix {
    std::vector<T> m_data;
public:
    const size_t m_nrow, m_ncol;
    Matrix(size_t nrow, size_t ncol):m_data(nrow*ncol, T(0)), m_nrow(nrow), m_ncol(ncol){}
    Matrix(size_t n): Matrix(n, n){}

    T& operator()(const size_t &irow, const size_t &icol) {
        ASSERT(irow >= 0 && irow < m_nrow);
        ASSERT(icol >= 0 && icol < m_ncol);
        return m_data[icol*m_nrow+irow];
    }

    bool is_square() const{
        return m_ncol==m_nrow;
    }

    EigenSolver<T> diagonalize() const{
        return EigenSolver<T>(*this);
    }

    void print()  {
        std::cout << std::setprecision(10);
        for (size_t irow=0ul; irow<m_nrow; ++irow){
            for (size_t icol=0ul; icol<m_ncol; ++icol){
                std::cout << operator()(irow, icol) << "  ";
            }
            std::cout << std::endl;
        }
    }

};

template<typename T>
class ColumnVector : public Matrix<T> {
public:
    ColumnVector(size_t nrow):Matrix<T>(nrow, 1){}
    T& operator()(const size_t &irow){
        return Matrix<T>::operator()(irow, 0);
    }
};

#endif //M7_MATRIX_H
