//
// Created by Robert J. Anderson on 09/01/2022.
//

#ifndef M7_ARNOLDISOLVER_H
#define M7_ARNOLDISOLVER_H

#include <M7_lib/linalg/DistMvProd.h>
#include <M7_lib/util/SmartPtr.h>

#include <arpackf.h>
#include <arrssym.h>
#include <arlnsmat.h>
#include <arrsnsym.h>
#include <arrscomp.h>

/**
 * put all type-independent operations in this un-templated base class
 */
struct ArnoldiProblemBase {
    const uint_t m_nroot;
    ArnoldiProblemBase(uint_t nroot): m_nroot(nroot){}
    /**
     * @return
     *  true if the obtained Krylov basis is sufficiently complete
     */
    virtual bool basis_found() = 0;

    /**
     * advance to the next step of the Arnoldi iteration
     */
    virtual void take_step() = 0;

    /**
     * @return
     *  true if another call to the Mv product provider is necessary
     */
    virtual bool do_another_mv_call() = 0;

    /**
     * after the Arnoldi iteration procedure has converged, prepare the eigenvalues
     */
    virtual void find_eigenvalues() = 0;

    /**
     * after the Arnoldi iteration procedure has converged, prepare the eigenvalues and the Ritz vectors
     */
    virtual void find_eigenvectors() = 0;

protected:
    bool solve_base(const std::function<void()> &product_fn, bool dist);
};

template<typename T>
struct ArnoldiProblemWithProduct : ArnoldiProblemBase {

    ArnoldiProblemWithProduct(uint_t nroot): ArnoldiProblemBase(nroot){}

    virtual void setup(uint_t nrow, bool parallel) = 0;

    virtual bool solve(dist_mv_prod::Base<T> &mv_prod) = 0;

    virtual bool solve(sparse::dynamic::Matrix<T> &sparse_mat) = 0;

    virtual T real_eigenvalue(uint_t i) = 0;

    virtual std::complex<T> complex_eigenvalue(uint_t i) = 0;

protected:
    /**
     * compute all required products in distributed mode
     * @param mv_prod
     *  MPI-distributed matrix-vector product
     */
    virtual void product(dist_mv_prod::Base<T> &mv_prod) = 0;

    /**
     * compute all required products in serial mode
     * @param sparse_mat
     *  sparse representation of the whole square matrix to be diagonalized on this MPI rank
     */
    virtual void product(sparse::dynamic::Matrix<T> &sparse_mat) = 0;
};


/**
 * Real, symmetric eigenvalue problem. It is the user's responsibility to ensure the operator is actually real symmetric
 * @tparam T
 *  floating point type of Krylov vectors
 */
template<typename T>
struct ArnoldiProblemSym : ArnoldiProblemWithProduct<T> {
    using ArnoldiProblemBase::m_nroot;
    using ArnoldiProblemBase::solve_base;

    std::unique_ptr<ARrcSymStdEig<T>> m_solver;

    ArnoldiProblemSym(uint_t nroot, T sigma=0.0) : ArnoldiProblemWithProduct<T>(nroot){
        dtype::unused(sigma);
    }

private:
    bool basis_found() override { return m_solver->ArnoldiBasisFound(); }

    void take_step() override { m_solver->TakeStep(); }

    bool do_another_mv_call() override {
        return (m_solver->GetIdo() == 1) || (m_solver->GetIdo() == -1);
    }

    void find_eigenvalues() override { m_solver->FindEigenvalues(); }

    void find_eigenvectors() override { m_solver->FindEigenvectors(); }

private:

    void setup(uint_t nrow, bool dist) override {
        if (mpi::i_am_root() || !dist) m_solver = smart_ptr::make_unique<ARrcSymStdEig<T>>(nrow, m_nroot);
    }

    void product(dist_mv_prod::Base<T> &mv_prod) override {
        mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
    }

    void product(sparse::dynamic::Matrix<T> &sparse_mat) override {
        sparse_mat.multiply(get_vector(), put_vector(), sparse_mat.nrow());
    }

    const T* get_vector() const {
        return m_solver ? m_solver->GetVector(): nullptr;
    }

    T* put_vector() const {
        return m_solver ? m_solver->PutVector(): nullptr;
    }

public:
    bool solve(dist_mv_prod::Base<T> &mv_prod) override {
        setup(mv_prod.m_nrow, true);
        auto prod_fn = [&]() {
            mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
        };
        return solve_base(prod_fn, true);
    }

    bool solve(sparse::dynamic::Matrix<T> &sparse_mat) override {
        setup(sparse_mat.nrow(), false);
        auto prod_fn = [&]() { sparse_mat.multiply(m_solver->GetVector(), m_solver->PutVector(), sparse_mat.nrow()); };
        return solve_base(prod_fn, false);
    }

    T real_eigenvalue(uint_t i) override {
        return m_solver->Eigenvalue(i);
    }

    std::complex<T> complex_eigenvalue(uint_t i) override {
        return {real_eigenvalue(i), 0.0};
    }
};

/**
 * Non-real, or non-symmetric eigenvalue problem.
 * @tparam T
 *  floating point type of Krylov vectors
 */
template<typename T>
struct ArnoldiProblemNonSym : ArnoldiProblemWithProduct<T> {
    using ArnoldiProblemBase::m_nroot;
    using ArnoldiProblemBase::solve_base;

    std::unique_ptr<ARrcNonSymStdEig<T>> m_solver;

    ArnoldiProblemNonSym(uint_t nroot, T sigma=0.0) : ArnoldiProblemWithProduct<T>(nroot){
        dtype::unused(sigma);
    }

private:
    bool basis_found() override { return m_solver->ArnoldiBasisFound(); }

    void take_step() override { m_solver->TakeStep(); }

    bool do_another_mv_call() override {
        return (m_solver->GetIdo() == 1) || (m_solver->GetIdo() == -1);
    }

    void find_eigenvalues() override { m_solver->FindEigenvalues(); }

    void find_eigenvectors() override { m_solver->FindEigenvectors(); }

    void setup(uint_t nrow, bool dist) override {
        if (mpi::i_am_root() || !dist) m_solver = smart_ptr::make_unique<ARrcNonSymStdEig<T>>(nrow, m_nroot);
    }

    void product(dist_mv_prod::Base<T> &mv_prod) override {
        mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
    }

    void product(sparse::dynamic::Matrix<T> &sparse_mat) override {
        sparse_mat.multiply(get_vector(), put_vector(), sparse_mat.nrow());
    }

    const T* get_vector() const {
        return m_solver ? m_solver->GetVector(): nullptr;
    }

    T* put_vector() const {
        return m_solver ? m_solver->PutVector(): nullptr;
    }

public:
    bool solve(dist_mv_prod::Base<T> &mv_prod) override {
        setup(mv_prod.m_nrow, true);
        auto prod_fn = [&]() {
            mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
        };
        return solve_base(prod_fn, true);
    }

    bool solve(sparse::dynamic::Matrix<T> &sparse_mat) override {
        setup(sparse_mat.nrow(), false);
        auto prod_fn = [&]() { sparse_mat.multiply(m_solver->GetVector(), m_solver->PutVector(), sparse_mat.nrow()); };
        return solve_base(prod_fn, false);
    }

    T real_eigenvalue(uint_t i) override {
        auto z = complex_eigenvalue(i);
        if (!fptol::numeric_zero(z.imag()))
            logging::warn("taking real part of eigenvalue with non-zero imaginary part");
        return arith::real(complex_eigenvalue(i));
    }

    std::complex<T> complex_eigenvalue(uint_t i) override {
        return m_solver->Eigenvalue(i);
    }
};


#endif //M7_ARNOLDISOLVER_H
