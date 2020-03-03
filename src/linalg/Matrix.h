//
// Created by Robert John Anderson on 2020-01-24.
//

#ifndef M7_MATRIX_H
#define M7_MATRIX_H

#include <iostream>
#include <vector>
#include <memory>
#include "../consts.h"
#include <Eigen/Eigenvalues>



template<typename T, bool herm>
struct eigensolver_t {
    typedef typename std::conditional<herm,
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic >>,
            Eigen::EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic >>>::type type;
};


template<typename T, bool herm>
struct eigensolver_t<std::complex<T>, herm> {
    typedef typename std::conditional<herm,
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic >>,
            Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic >>>::type type;
};

template <typename T, bool herm>
class Matrix {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic > M_t;
    typedef typename eigensolver_t<T, herm>::type ES_t;
public:
    const size_t m_nrow, m_ncol;
private:
    M_t m_data;
public:
    Matrix(const size_t &nrow, const size_t &ncol):
    m_nrow(nrow), m_ncol(ncol), m_data(m_nrow, m_ncol){}

    Matrix(const size_t &nrow): Matrix(nrow, nrow){}


    T get(const size_t &irow, const size_t &icol){
        return m_data(irow, icol);
    }
    void set(const size_t &irow, const size_t &icol, const T &value){
        m_data(irow, icol) = value;
        if (herm && irow!=icol) {
            assert(consts::float_is_zero(get(icol, irow)) ||
                    consts::floats_nearly_equal(get(icol, irow), consts::conj(value)));
            m_data(icol, irow) = consts::conj(value);
        }
    }
    bool is_square(){
        return m_nrow==m_ncol;
    }
    std::unique_ptr<ES_t> diagonalise(){
        assert(is_square());
        auto es = std::make_unique<ES_t>();
        es.get()->compute(m_data);
        return es;
    }
};

#endif //M7_MATRIX_H
