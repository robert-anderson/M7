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
 * options to pass to the ARPACK solver
 */
struct ArnoldiOptions {
    /**
     * number of eigenpairs be computed
     */
    uint_t m_nroot = 2ul;
    /**
     * number of arnoldi vectors to be generated at each iteration
     */
    uint_t m_narnoldi_vector = 0ul;
    /**
     * maximum iteration number
     */
    uint_t m_niter_max = 0ul;
    /**
     * ritz vector tolerance determining convergence criterion
     */
    double m_ritz_tol = 0.0;
};

/**
 * put all type-independent operations in this un-templated base class
 */
struct ArnoldiProblemBase {

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
    virtual bool find_eigenvectors() = 0;

protected:
    bool solve_base(const std::function<void()> &product_fn, bool dist);

    template<class arfloat_t, class ar_t>
    void begin_log(ARrcStdEig<arfloat_t, ar_t>* solver, bool sym) {
        REQUIRE_TRUE(solver, "solver object must be allocated");
        logging::info("Finding extremal eigenpairs with {}symmetric Arnoldi method...", sym ? "" : "non-");
        logging::info("Total eigenproblem dimension: {}", solver->GetN());
        logging::info("Number of eigenpairs to converge: {}", solver->GetNev());
        logging::info("Ritz vector convergence threshold: {}", solver->GetTol());
        logging::info("Maximum number of Arnoldi iterations: {}", solver->GetMaxit());
        logging::info("Number of Arnoldi vectors generated at each iteration: {}", solver->GetNcv());
    }

    template<class arfloat_t, class ar_t>
    void end_log(ARrcStdEig<arfloat_t, ar_t>* solver, bool success) {
        REQUIRE_TRUE(solver, "solver object must be allocated");
        const uint_t niter = solver->GetIter();
        logging::info("Arnoldi {} after {} iteration{}",
                      success ? "converged" : "failed to converge",
                      niter, string::plural(niter));
    }
};

template<typename T>
struct ArnoldiProblemWithProduct : ArnoldiProblemBase {


    virtual void setup(uint_t nrow, ArnoldiOptions opts, bool parallel) = 0;

    /**
     * distributed MV product version of the Arnoldi solver
     * @param mv_prod
     *  object providing access to the MV product distributed over all processes
     * @param opts
     *  configuring options of the ARPACK interface
     * @return
     *  true if all eigenpairs converged
     */
    virtual bool solve(dist_mv_prod::Base<T> &mv_prod, ArnoldiOptions opts) = 0;

    /**
     * non-distributed MV product version of the Arnoldi solver
     * @param sparse_mat
     *  sparse matrix object held only on this MPI rank
     * @param opts
     *  configuring options of the ARPACK interface
     * @return
     *  true if all eigenpairs converged
     */
    virtual bool solve(sparse::dynamic::Matrix<T> &sparse_mat, ArnoldiOptions opts) = 0;

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
    using ArnoldiProblemBase::solve_base;

    std::unique_ptr<ARrcSymStdEig<T>> m_solver;

private:
    bool basis_found() override { return m_solver->ArnoldiBasisFound(); }

    void take_step() override { m_solver->TakeStep(); }

    bool do_another_mv_call() override {
        return (m_solver->GetIdo() == 1) || (m_solver->GetIdo() == -1);
    }

    void find_eigenvalues() override { m_solver->FindEigenvalues(); }

    bool find_eigenvectors() override {
        m_solver->FindEigenvectors();
        return m_solver->EigenvectorsFound();
    }

private:

    void setup(uint_t nrow, ArnoldiOptions opts, bool dist) override {
        if (mpi::i_am_root() || !dist) {
            m_solver = smart_ptr::make_unique<ARrcSymStdEig<T>>(
                    nrow, opts.m_nroot, "LM", opts.m_narnoldi_vector, opts.m_ritz_tol, opts.m_niter_max);
            ArnoldiProblemBase::begin_log(m_solver.get(), true);
        }
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
    bool solve(dist_mv_prod::Base<T> &mv_prod, ArnoldiOptions opts) override {
        setup(mv_prod.m_nrow, opts, true);
        auto prod_fn = [&]() {
            mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
        };
        const auto success = solve_base(prod_fn, true);
        if (mpi::i_am_root()) ArnoldiProblemBase::end_log(m_solver.get(), success);
        return success;
    }

    bool solve(sparse::dynamic::Matrix<T> &sparse_mat, ArnoldiOptions opts) override {
        setup(sparse_mat.nrow(), opts, false);
        auto prod_fn = [&]() { sparse_mat.multiply(m_solver->GetVector(), m_solver->PutVector(), sparse_mat.nrow()); };
        return solve_base(prod_fn, false);
    }

    T real_eigenvalue(uint_t i) override {
        i = (m_solver->GetNev()-i)-1;
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
    using ArnoldiProblemBase::solve_base;

    std::unique_ptr<ARrcNonSymStdEig<T>> m_solver;

private:
    bool basis_found() override { return m_solver->ArnoldiBasisFound(); }

    void take_step() override { m_solver->TakeStep(); }

    bool do_another_mv_call() override {
        return (m_solver->GetIdo() == 1) || (m_solver->GetIdo() == -1);
    }

    void find_eigenvalues() override { m_solver->FindEigenvalues(); }

    bool find_eigenvectors() override {
        m_solver->FindEigenvectors();
        return m_solver->EigenvectorsFound();
    }

    void setup(uint_t nrow, ArnoldiOptions opts, bool dist) override {
        if (mpi::i_am_root() || !dist) {
            m_solver = smart_ptr::make_unique<ARrcNonSymStdEig<T>>(
                    nrow, opts.m_nroot, "LM", opts.m_narnoldi_vector, opts.m_ritz_tol, opts.m_niter_max);
            ArnoldiProblemBase::begin_log(m_solver.get(), false);
        }
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
    bool solve(dist_mv_prod::Base<T> &mv_prod, ArnoldiOptions opts) override {
        setup(mv_prod.m_nrow, opts, true);
        auto prod_fn = [&]() {
            mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
        };
        const auto success = solve_base(prod_fn, true);
        if (mpi::i_am_root()) ArnoldiProblemBase::end_log(m_solver.get(), success);
        return success;
    }

    bool solve(sparse::dynamic::Matrix<T> &sparse_mat, ArnoldiOptions opts) override {
        setup(sparse_mat.nrow(), opts, false);
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
        /*
         * nonsym state ordering seems to be opposite to that of sym
         */
        //i = (m_solver->GetNev()-i)-1;
        return m_solver->Eigenvalue(i);
    }
};


#endif //M7_ARNOLDISOLVER_H
