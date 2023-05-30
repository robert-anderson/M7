//
// Created by anderson on 18/07/2022.
//

#ifndef M7_CIINITIALIZER_H
#define M7_CIINITIALIZER_H

#include <M7_lib/linalg/FciIters.h>
#include <M7_lib/arnoldi/ArnoldiSolver.h>

namespace ci_init {
    struct Options : ArnoldiOptions {
        /**
         * shift to add to the diagonal elements of the sparse subspace Hamiltonian
         */
        ham_comp_t m_diag_shift = 0.0;
        /**
         * enumerating loop kind: connections from each MBF or MBF pairs.
         */
        enum LoopKind {Conns, MbfPairs};
        LoopKind m_loop_kind = Conns;
    };

    struct Initializer {
        const Options m_opts;
        const bool m_is_hermitian;
        sparse::dynamic::Matrix<ham_t> m_sparse_ham;

        /**
         * mapped list of basis functions to aid in the setup of sparse H, and retain the physical meaning of its rows
         */
        typedef SingleFieldRow<field::Mbf> mbf_order_row_t;
        typedef buffered::MappedTable<mbf_order_row_t> mbf_order_table_t;
        mbf_order_table_t m_mbf_order_table;

        Initializer(const Hamiltonian& h, sys::Particles particles, Options opts = {});

        explicit Initializer(const Hamiltonian& h, Options opts = {});

    private:

        /**
         * build Hamiltonian in subspace by looping over mbfs and then by connections (recommended for large spaces)
         */
        void build_ham_conns(const Hamiltonian &h, ham_comp_t diag_shift);

        /**
         * build Hamiltonian in subspace by looping over pairs of mbfs (recommended for small spaces)
         */
        void build_ham_mbfs(const Hamiltonian &h, ham_comp_t diag_shift);

        template<uint_t sym>
        ArnoldiSolver<ham_t> solve(tag::Int<sym>) {
            dist_mv_prod::Sparse<ham_t> dist(m_sparse_ham);
            ArnoldiSolver<ham_t> solver(dist, m_opts, tag::Int<sym>());
            /*
             * once the ARPACK procedure is complete, the eigenvalues must be adjusted to undo the diagonal shift
             */
            solver.shift_evals(-m_opts.m_diag_shift);
            return solver;
        }

    public:

        ArnoldiSolver<ham_t> solve() {
            return m_is_hermitian ? solve(ArnoldiSolverBase::c_sym) : solve(ArnoldiSolverBase::c_nonsym);
        }

        /**
         * in instances where retention of the MBF list and sparse Hamiltonian is not desired
         */
        static ArnoldiSolver<ham_t> solve(const Hamiltonian& h, sys::Particles particles, Options opts = {}) {
            return Initializer(h, particles, opts).solve();
        }

        static ArnoldiSolver<ham_t> solve(const Hamiltonian& h, Options opts = {}) {
            return Initializer(h, h.default_particles(), opts).solve();
        }
    };
}

#endif //M7_CIINITIALIZER_H
