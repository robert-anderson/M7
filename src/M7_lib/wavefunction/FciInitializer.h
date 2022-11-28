//
// Created by anderson on 18/07/2022.
//

#ifndef M7_FCIINITIALIZER_H
#define M7_FCIINITIALIZER_H

#include <M7_lib/linalg/FciIters.h>
#include <M7_lib/arnoldi/ArnoldiSolver.h>

struct FciInitOptions : ArnoldiOptions {
    /**
     * shift to add to the diagonal elements of the sparse subspace Hamiltonian
     */
    ham_comp_t m_diag_shift = 0.0;
};

struct FciInitializer {
    const FciInitOptions m_opts;
    const bool m_is_hermitian;
    sparse::dynamic::Matrix<ham_t> m_sparse_ham;

    struct MbfOrderRow : Row {
        field::Mbf m_mbf;
        MbfOrderRow(sys::Basis basis): m_mbf(this, basis, "key"){}

        field::Mbf &key_field(){return m_mbf;}
    };
    /**
     * mapped list of basis functions to aid in the setup of sparse H, and retain the physical meaning of its rows
     */
    buffered::MappedTable<MbfOrderRow> m_mbf_order_table;

    explicit FciInitializer(const Hamiltonian& h, FciInitOptions opts={});

private:
    template<uint_t sym>
    ArnoldiSolver<ham_t> solve(tag::Int<sym>){
        dist_mv_prod::Sparse<ham_t> dist(m_sparse_ham);
        ArnoldiSolver<ham_t> solver(dist, m_opts, tag::Int<sym>());
        /*
         * once the ARPACK procedure is complete, the eigenvalues must be adjusted to undo the diagonal shift
         */
        solver.shift_evals(-m_opts.m_diag_shift);
        return solver;
    }
public:

    ArnoldiSolver<ham_t> solve(){
        return m_is_hermitian ? solve(ArnoldiSolverBase::c_sym) : solve(ArnoldiSolverBase::c_nonsym);
    }

    /**
     * in instances where retention of the MBF list and sparse Hamiltonian is not desired
     */
    static ArnoldiSolver<ham_t> solve(const Hamiltonian& h, FciInitOptions opts={}) {
        return FciInitializer(h, opts).solve();
    }
};


#endif //M7_FCIINITIALIZER_H
