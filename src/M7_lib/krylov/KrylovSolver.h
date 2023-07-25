//
// Created by rja on 25/07/23.
//

#ifndef M7_KRYLOVSOLVER_H
#define M7_KRYLOVSOLVER_H

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/Arith.h>
#include <algorithm>
#include <numeric>

struct KrylovOptions {
    /**
     * maximum iteration number
     */
    uint_t m_niter_max = 0ul;
    /**
     * ritz vector tolerance determining convergence criterion
     */
    double m_ritz_tol = 1e-9;
};


template<typename kry_t>
struct KrylovSolver {
    typedef arith::comp_t<kry_t> comp_t;
    /*
     * eigenvectors are by definition kry_t, but eigenvalues can be real or complex depending on symmetry, so their
     * real/imag parts are stored separately
     */

    const uint_t m_nroot;
    const uint_t m_nelement_evec;

    KrylovSolver(uint_t nroot, uint_t nelement_evec): m_nroot(nroot), m_nelement_evec(nelement_evec){}
protected:
    uintv_t m_root_ordering;
    v_t<comp_t> m_real_evals;
    v_t<comp_t> m_imag_evals;
    v_t<kry_t> m_evecs;

    void set_results(const comp_t* real_evals, const comp_t* imag_evals, const kry_t* raw_evecs){
        if (!real_evals) return;

        m_real_evals = {real_evals, real_evals + m_nroot};
        if (imag_evals) m_imag_evals = {imag_evals, imag_evals + m_nroot};
        else m_imag_evals.assign(m_nroot, 0.0);
        m_evecs = {raw_evecs, raw_evecs+(m_nroot * m_nelement_evec)};
        m_root_ordering.resize(m_nroot);
        std::iota(m_root_ordering.begin(), m_root_ordering.end(), 0);
        // sort with largest-magnitude eval first
        std::sort(m_root_ordering.begin(), m_root_ordering.end(), [&](uint_t i, uint_t j) {
            std::complex<comp_t> zi = {m_real_evals[i], m_imag_evals[i]};
            std::complex<comp_t> zj = {m_real_evals[j], m_imag_evals[j]};
            return std::abs(zi) > std::abs(zj);
        });
    }

    /**
     * send the eigenvalues to each process
     */
    void bcast(uint_t irank=0ul) {
        mpi::bcast(m_root_ordering, irank);
        mpi::bcast(m_real_evals, irank);
        mpi::bcast(m_imag_evals, irank);
    }

public:

    void get_eval(uint_t iroot, comp_t& eval) const {
        REQUIRE_NEAR_ZERO(m_imag_evals[m_root_ordering[iroot]], "non-zero imaginary part");
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
        for (size_t iroot = 0ul; iroot < m_nroot; ++iroot) {
            evals.push_back({});
            get_eval(m_root_ordering[iroot], evals.back());
        }
    }

public:

    void shift_evals(comp_t shift) {
        if (m_real_evals.empty()) return;
        for (auto& it: m_real_evals) it+=shift;
    }

    void get_evals(v_t<comp_t>& evals) const { get_evals<true>(evals, 0); }
    void get_evals(v_t<std::complex<comp_t>>& evals) const { get_evals<false>(evals, 0); }

    const kry_t* get_evec(uint_t iroot) const {
        if (m_evecs.empty()) return nullptr;
        REQUIRE_LT(iroot, m_nroot, "root index OOB");
        return m_evecs.data()+(m_root_ordering[iroot] * m_nelement_evec);
    }

    v_t<const kry_t*> get_evecs() const {
        if (m_evecs.empty()) return {};
        v_t<const kry_t*> evecs;
        for (uint_t iroot = 0ul; iroot < m_nroot; ++iroot) evecs.push_back(get_evec(iroot));
        return evecs;
    }

    bool i_have_evecs() const {
        return !m_evecs.empty();
    }
};


#endif //M7_KRYLOVSOLVER_H
