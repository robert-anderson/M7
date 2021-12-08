//
// Created by Robert John Anderson on 2020-01-24.
//

#include "Matrix.h"

template<>
void Matrix<double>::multiply(const std::vector<double> &in, std::vector<double> &out) {
    char trans = 'N';
    const double alpha = 1.0;
    const double beta = 1.0;
    int nrow = m_nrow;
    int ncol = m_ncol;
    int incx = 1;
    int incy = 1;
    dgemv_(&trans, &ncol, &nrow, &alpha, m_buffer.data(), &nrow, in.data(), &incx, &beta, out.data(), &incy);
}

template<>
void Matrix<std::complex<double>>::multiply(const std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out) {
    char trans = 'N';
    const std::complex<double> alpha = 1.0;
    const std::complex<double> beta = 1.0;
    int nrow = m_nrow;
    int ncol = m_ncol;
    int incx = 1;
    int incy = 1;
    zgemv_(&trans, &ncol, &nrow, &alpha, m_buffer.data(), &nrow, in.data(), &incx, &beta, out.data(), &incy);
}