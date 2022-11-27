//
// Created by Robert J. Anderson on 09/01/2022.
//

#ifndef M7_ARNOLDISOLVER_H
#define M7_ARNOLDISOLVER_H

#include <M7_lib/linalg/DistMvProd.h>
#include <M7_lib/util/Pointer.h>

#include <arpackf.h>
#include <arrssym.h>
#include <arlnsmat.h>
#include <arrsnsym.h>
#include <arrscomp.h>
#include <numeric>
#include <M7_lib/util/Sort.h>


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
    double m_ritz_tol = 1e-9;
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
struct ArnoldiResults {
    static_assert(!dtype::is_complex<T>(), "Results are held in component type");
protected:
    v_t<T> m_real_evals;
    v_t<T> m_imag_evals;
    /**
     * vector of pointers to the real part of the first element of each eigenvector
     */
    v_t<const T*> m_evecs;
    /**
     * number of elements in each evec
     */
    uint_t m_nelement_evec = 0ul;
    /**
     * true if the eigenvectors are complex-valued (vector is in the real_0, imag_0, real_1, imag_1, ... pattern)
     */
    bool m_complex_evecs = false;
public:
    ArnoldiResults(uint_t nroot, uint_t nelement_evec, const T* real_evals, const T* imag_evals,
                   const T* raw_evecs, bool complex_evecs):
        m_nelement_evec(nelement_evec), m_complex_evecs(complex_evecs) {
        if (!real_evals) return;

        m_real_evals = v_t<T>(real_evals, real_evals + nroot);
        if (imag_evals) m_imag_evals = v_t<T>(imag_evals, imag_evals + nroot);
        else m_imag_evals.assign(nroot, 0.0);
        for (uint_t i = 0; i < nroot; ++i)
            m_evecs.push_back(raw_evecs + i * m_nelement_evec * (m_complex_evecs ? 2 : 1));
        uintv_t ordering(nroot);
        std::iota(ordering.begin(), ordering.end(), 0);
        // sort with largest-magnitude eval first
        std::sort(ordering.begin(), ordering.end(), [&](uint_t i, uint_t j) {
            std::complex<T> zi = {m_real_evals[i], m_imag_evals[i]};
            std::complex<T> zj = {m_real_evals[j], m_imag_evals[j]};
            return std::abs(zi) > std::abs(zj);
        });
        m_real_evals = sort::reorder(m_real_evals, ordering);
        m_imag_evals = sort::reorder(m_imag_evals, ordering);
        m_evecs = sort::reorder(m_evecs, ordering);
    }

    ArnoldiResults(){}

    void get_eval(uint_t iroot, T& eval) const {
        REQUIRE_FALSE(m_imag_evals[iroot], "non-zero imaginary part");
        eval = m_real_evals[iroot];
    }

    void get_eval(uint_t iroot, std::complex<T>& eval) const {
        eval = {m_real_evals[iroot], m_imag_evals[iroot]};
    }

private:
    /*
     * final arg is a dummy to enable static dispatch in the arithmetic-resolving methods below
     */
    template<bool real>
    void get_evals(v_t<arith::num_t<T, real>>& evals, int) const {
        evals.clear();
        for (size_t iroot = 0ul; iroot < nroot(); ++iroot) {
            evals.push_back({});
            get_eval(iroot, evals.back());
        }
    }

public:
    void get_evals(v_t<T>& evals) const { get_evals<true>(evals, 0); }
    void get_evals(v_t<std::complex<T>>& evals) const { get_evals<false>(evals, 0); }

    uint_t nroot() const {
        return m_real_evals.size();
    }

    uint_t nelement_evec() const {
        return m_nelement_evec;
    }

private:
    template<bool real>
    void get_evec(uint_t iroot, const arith::num_t<T, real>* &evec, int) const {
        if (real) {
            REQUIRE_FALSE(m_complex_evecs, "cannot dereference complex eigenvector as a real vector");
        }
        else {
            REQUIRE_TRUE(m_complex_evecs, "cannot dereference real eigenvector as a complex vector");
        }
        if (m_evecs.empty()) {
            // this rank does not have access to the eigenvectors
            evec = nullptr;
            return;
        }
        REQUIRE_LT(iroot, nroot(), "root index OOB");
        evec = reinterpret_cast<const arith::num_t<T, real>*>(m_evecs[iroot]);
    }

public:
    void get_evec(uint_t iroot, const T* &evec) const { get_evec<true>(iroot, evec, 0);}
    void get_evec(uint_t iroot, const std::complex<T>* &evec) const { get_evec<false>(iroot, evec, 0);}

private:
    template<bool real>
    void get_evecs(v_t<const arith::num_t<T, real>*>& evecs, int) const {
        if (m_evecs.empty()) {
            // this rank does not have access to the eigenvectors
            evecs.clear();
            return;
        }
        evecs.clear();
        for (uint_t iroot = 0ul; iroot<nroot(); ++iroot) {
            evecs.push_back(nullptr);
            get_evec(iroot, evecs.back());
        }
    }

public:
    void get_evecs(v_t<const T*>& evecs) const { get_evecs<true>(evecs, 0);}
    void get_evecs(v_t<const std::complex<T>*>& evecs) const { get_evecs<false>(evecs, 0);}

    /**
     * send the eigenvalues to each process
     */
    void bcast(uint_t irank=0ul) {
        mpi::bcast(m_real_evals, irank);
        mpi::bcast(m_imag_evals, irank);
        mpi::bcast(m_nelement_evec, irank);
        mpi::bcast(m_complex_evecs, irank);
    }

    bool i_have_evecs() const {
        return !m_evecs.empty();
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

    virtual ArnoldiResults<arith::comp_t<T>> get_results() = 0;

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
        return m_solver->FindEigenvectors();
    }

private:

    void setup(uint_t nrow, ArnoldiOptions opts, bool dist) override {
        if (mpi::i_am_root() || !dist) {
            m_solver = ptr::smart::make_unique<ARrcSymStdEig<T>>(
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

    ArnoldiResults<arith::comp_t<T>> get_results() override {
        if (!m_solver) return {};
        return {uint_t(m_solver->GetNev()), uint_t(m_solver->GetN()),
                m_solver->RawEigenvalues(), nullptr, m_solver->RawEigenvectors(), false};
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
        return m_solver->FindEigenvectors();
    }

    void setup(uint_t nrow, ArnoldiOptions opts, bool dist) override {
        if (mpi::i_am_root() || !dist) {
            m_solver = ptr::smart::make_unique<ARrcNonSymStdEig<T>>(
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

    ArnoldiResults<arith::comp_t<T>> get_results() override {
        if (!m_solver) return {};
        return {uint_t(m_solver->GetNev()), uint_t(m_solver->GetN()),
                m_solver->RawEigenvalues(), m_solver->RawEigenvaluesImag(), m_solver->RawEigenvectors(), false};
    }
};





/**
 * put all type-independent operations in this un-templated base class
 */
struct ArnoldiSolverBase {
protected:
    const uint_t m_nroot;
    const uint_t m_nelement_evec;

public:

    /**
     * tags to statically specify symmetric or non-symmetric ARPACK algorithms
     */
    static constexpr tag::Int<1> c_sym = {};
    static constexpr tag::Int<0> c_nonsym = {};

    ArnoldiSolverBase(uint_t nroot, uint_t nelement_evec): m_nroot(nroot), m_nelement_evec(nelement_evec){}
    template<typename comp_t, bool real, bool sym> struct SolverSelector {};

    template<typename comp_t> struct SolverSelector<comp_t, true, true>{ typedef ARrcSymStdEig<comp_t> type;};
    template<typename comp_t> struct SolverSelector<comp_t, true, false>{ typedef ARrcNonSymStdEig<comp_t> type;};
    template<typename comp_t> struct SolverSelector<comp_t, false, true>{ typedef ARrcCompStdEig<comp_t> type;};
    template<typename comp_t> struct SolverSelector<comp_t, false, false>{ typedef ARrcCompStdEig<comp_t> type;};

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
     * after the Arnoldi iteration procedure has converged, prepare the eigenvalues and the Ritz vectors
     */
    virtual bool find_eigenvectors() = 0;

protected:
    bool solve(const std::function<void()> &product_fn, bool dist);

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

template<typename kry_t>
struct ArnoldiSolver : ArnoldiSolverBase {
    typedef arith::comp_t<kry_t> comp_t;
    /*
     * eigenvectors are by definition kry_t, but eigenvalues can be real or complex depending on symmetry, so their
     * real/imag parts are stored separately
     */
protected:
    uintv_t m_root_ordering;
    v_t<comp_t> m_real_evals;
    v_t<comp_t> m_imag_evals;
    v_t<kry_t> m_evecs;
    ARrcStdEig<comp_t, kry_t>* m_ar_base = nullptr;

private:
    void set_results(const comp_t* real_evals, const comp_t* imag_evals, const kry_t* raw_evecs){
        if (!real_evals) return;

        m_real_evals = {real_evals, real_evals + m_nroot};
        if (imag_evals) m_imag_evals = {imag_evals, imag_evals + m_nroot};
        else m_imag_evals.assign(m_nroot, 0.0);
        m_evecs = {raw_evecs, raw_evecs+(m_nroot*m_nelement_evec)};
        m_root_ordering.resize(m_nroot);
        std::iota(m_root_ordering.begin(), m_root_ordering.end(), 0);
        // sort with largest-magnitude eval first
        std::sort(m_root_ordering.begin(), m_root_ordering.end(), [&](uint_t i, uint_t j) {
            std::complex<comp_t> zi = {m_real_evals[i], m_imag_evals[i]};
            std::complex<comp_t> zj = {m_real_evals[j], m_imag_evals[j]};
            return std::abs(zi) > std::abs(zj);
        });
    }

    void set_results(ARrcSymStdEig<comp_t>* ar) {
        if (!ar) return;
        set_results(ar->RawEigenvalues(), nullptr, ar->RawEigenvectors());
    }

    void set_results(ARrcNonSymStdEig<comp_t>* ar) {
        if (!ar) return;
        set_results(ar->RawEigenvalues(), ar->RawEigenvaluesImag(), ar->RawEigenvectors());
    }

    void set_results(ARrcCompStdEig<comp_t>* ar) {
        if (!ar) return;
        v_t<comp_t> evals_re(m_nroot);
        v_t<comp_t> evals_im(m_nroot);
        for (uint_t iroot=0ul; iroot<m_nroot; ++iroot){
            const auto eval = ar->Eigenvalue(iroot);
            evals_re[iroot] = eval.real();
            evals_im[iroot] = eval.imag();
        }
        set_results(evals_re.data(), evals_im.data(), ar->RawEigenvectors());
    }

    /**
     * send the eigenvalues to each process
     */
    void bcast(uint_t irank=0ul) {
        mpi::bcast(m_root_ordering, irank);
        mpi::bcast(m_real_evals, irank);
        mpi::bcast(m_imag_evals, irank);
    }

    ArnoldiSolver(uint_t nroot, uint_t nelement_evec): ArnoldiSolverBase(nroot, nelement_evec){}

    template<uint_t sym>
    ArnoldiSolver(std::function<void()> prod_fn, bool dist, uint_t nelement_evec, ArnoldiOptions opts, tag::Int<sym>):
            ArnoldiSolver(opts.m_nroot, nelement_evec) {
        constexpr bool real = !dtype::is_complex<kry_t>();
        typedef typename ArnoldiSolverBase::SolverSelector<comp_t, real, sym>::type ar_t;
        std::unique_ptr<ar_t> m_ar;
        if (mpi::i_am_root() || !dist) {
            m_ar = ptr::smart::make_unique<ar_t>(
                    m_nelement_evec, opts.m_nroot, "LM", opts.m_narnoldi_vector, opts.m_ritz_tol, opts.m_niter_max);
            //ArnoldiProblemBase::begin_log(m_solver.get(), false);
            m_ar_base = m_ar.get();
        }

        const auto success = ArnoldiSolverBase::solve(prod_fn, dist);
        //if (mpi::i_am_root()) ArnoldiProblemBase::end_log(m_solver.get(), success);
        if (success) set_results(m_ar.get());
        bcast();
    }

    /*
     * get and put vectors only exist on ranks with a solver instance
     */
    const kry_t* get_vector() {
        return m_ar_base ? m_ar_base->GetVector() : nullptr;
    }
    kry_t* put_vector() {
        return m_ar_base ? m_ar_base->PutVector() : nullptr;
    }

public:
    template<uint_t sym>
    ArnoldiSolver(sparse::dynamic::Matrix<kry_t> &sparse_mat, uint_t nelement_evec, ArnoldiOptions opts, tag::Int<sym>):
            ArnoldiSolver([&](){
                sparse_mat.multiply(get_vector(), put_vector(), nelement_evec);
            }, false, nelement_evec, opts, tag::Int<sym>()) {}

    template<uint_t sym>
    ArnoldiSolver(dist_mv_prod::Base<kry_t> &mv_prod, ArnoldiOptions opts, tag::Int<sym>):
            ArnoldiSolver([&](){
                mv_prod.parallel_multiply(get_vector(), mv_prod.m_nrow, put_vector());
            }, true, mv_prod.m_nrow, opts, tag::Int<sym>()) {}

    bool basis_found() override { return m_ar_base->ArnoldiBasisFound(); }

    void take_step() override { m_ar_base->TakeStep(); }

    bool do_another_mv_call() override {
        return (m_ar_base->GetIdo() == 1) || (m_ar_base->GetIdo() == -1);
    }

    bool find_eigenvectors() override {
        return m_ar_base->FindEigenvectors();
    }

    void get_eval(uint_t iroot, comp_t& eval) const {
        REQUIRE_FALSE(m_imag_evals[m_root_ordering[iroot]], "non-zero imaginary part");
        eval = m_real_evals[m_root_ordering[iroot]];
    }

    void get_eval(uint_t iroot, std::complex<comp_t>& eval) const {
        eval = {m_real_evals[m_root_ordering[iroot]], m_imag_evals[m_root_ordering[iroot]]};
    }

private:
    /*
     * final arg is a dummy to enable static dispatch in the arithmetic-resolving methods below
     */
    template<bool real>
    void get_evals(v_t<arith::num_t<comp_t, real>>& evals, int) const {
        evals.clear();
        for (size_t iroot = 0ul; iroot < nroot(); ++iroot) {
            evals.push_back({});
            get_eval(m_root_ordering[iroot], evals.back());
        }
    }

public:
    void get_evals(v_t<comp_t>& evals) const { get_evals<true>(evals, 0); }
    void get_evals(v_t<std::complex<comp_t>>& evals) const { get_evals<false>(evals, 0); }

    uint_t nroot() const {
        return m_real_evals.size();
    }

    uint_t nelement_evec() const {
        return m_nelement_evec;
    }

    void get_evec(uint_t iroot, const kry_t* &evec) const {
        if (m_evecs.empty()) {
            // this rank does not have access to the eigenvectors
            evec = nullptr;
            return;
        }
        REQUIRE_LT(iroot, nroot(), "root index OOB");
        evec = m_evecs.data()+(m_root_ordering[iroot]*m_nelement_evec);
    }

    template<bool real>
    void get_evecs(v_t<const kry_t*>& evecs) const {
        if (m_evecs.empty()) {
            // this rank does not have access to the eigenvectors
            evecs.clear();
            return;
        }
        evecs.clear();
        for (uint_t iroot = 0ul; iroot<nroot(); ++iroot) {
            evecs.push_back(nullptr);
            get_evec(iroot, evecs.back());
        }
    }

    bool i_have_evecs() const {
        return !m_evecs.empty();
    }
};

#endif //M7_ARNOLDISOLVER_H
