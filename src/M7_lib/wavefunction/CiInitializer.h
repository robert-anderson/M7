//
// Created by anderson on 18/07/2022.
//

#ifndef M7_CIINITIALIZER_H
#define M7_CIINITIALIZER_H

#include <M7_lib/linalg/FciIters.h>
#include <M7_lib/arnoldi/ArnoldiSolver.h>

namespace ci_init {
    struct Options {
        /**
         * number of eigenpairs be computed
         */
        uint_t m_nroot = 1ul;
        /**
         * shift to add to the diagonal elements of the sparse subspace Hamiltonian
         */
        ham_comp_t m_diag_shift = 0.0;
        /**
         * enumerating loop kind: connections from each MBF or MBF pairs.
         */
        enum LoopKind {Conns, MbfPairs};
        LoopKind m_loop_kind = Conns;

        enum SolverKind {NoSolve, Davidson, Arnoldi};
        SolverKind m_solver_kind = NoSolve;
    };

    struct Subspace {
        /**
         * mapped list of basis functions to aid in the setup of sparse H, and retain the physical meaning of its rows
         */
        mbf::table_t m_mbf_order_table;
        const Hamiltonian* m_h;
        Subspace(const Hamiltonian* h);
    };

    struct FciSubspace : Subspace {
        FciSubspace(const Hamiltonian* h, sys::Particles particles);
        explicit FciSubspace(const Hamiltonian* h): FciSubspace(h, h->default_particles()){}
    };

    struct RefConnSubspace : Subspace {
        RefConnSubspace(const Hamiltonian* h, const Mbf& ref);
    };

    struct Initializer {
        const Options m_opts;
        const bool m_is_hermitian;
        sparse::dynamic::Matrix<ham_t> m_sparse_ham;

        Initializer(const Subspace& subspace, Options opts = {});

    private:

        /**
         * build Hamiltonian in subspace by looping over mbfs and then by connections (recommended for large spaces)
         */
        void build_ham_conns(const Subspace& subspace, ham_comp_t diag_shift);

        /**
         * build Hamiltonian in subspace by looping over pairs of mbfs (recommended for small spaces)
         */
        void build_ham_mbfs(const Subspace& subspace, ham_comp_t diag_shift);

#ifdef ENABLE_ARPACK
        template<uint_t sym>
        ArnoldiSolver<ham_t> solve_arnoldi(tag::Int<sym>) {
            dist_mv_prod::Sparse<ham_t> dist(m_sparse_ham);
            ArnoldiSolver<ham_t> solver(dist, m_opts, tag::Int<sym>());
            /*
             * once the ARPACK procedure is complete, the eigenvalues must be adjusted to undo the diagonal shift
             */
            solver.shift_evals(-m_opts.m_diag_shift);
            return solver;
        }
#endif

    public:

#ifdef ENABLE_ARPACK
        ArnoldiSolver<ham_t> solve_arnoldi() {
            return m_is_hermitian ? solve(ArnoldiSolverBase::c_sym) : solve(ArnoldiSolverBase::c_nonsym);
        }

        /**
         * in instances where retention of the sparse Hamiltonian is not desired
         */
        static ArnoldiSolver<ham_t> solve(const Subspace& subspace, Options opts = {}) {
            return Initializer(subspace, opts).solve_arnoldi();
        }
#endif
    };
}

#endif //M7_CIINITIALIZER_H
