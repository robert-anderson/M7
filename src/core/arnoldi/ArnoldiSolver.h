//
// Created by rja on 09/01/2022.
//

#ifndef M7_ARNOLDISOLVER_H
#define M7_ARNOLDISOLVER_H

#include "src/core/linalg/DistMvProd.h"
#include "arpackpp/include/arpackf.h"
#include "arpackpp/include/arrssym.h"
#include "arpackpp/include/arlnsmat.h"
#include "arpackpp/include/arrsnsym.h"
#include "arpackpp/include/arrscomp.h"

template<typename T, typename U>
struct ArnoldiSolver {
//    std::unique_ptr<ARrcSymStdEig<double>> m_arpack_solver_sym;
    std::unique_ptr<ARrcNonSymStdEig<double>> m_arpack_solver_nonsym;
//    std::unique_ptr<ARrcSymStdEig<double>> m_arpack_solver_complex;


    bool solve(dist_mv_prod::Base<T>& mat, size_t nroot) {
        m_arpack_solver_nonsym = std::unique_ptr<ARrcNonSymStdEig<double>>(new ARrcNonSymStdEig<double> (mat.m_nrow, nroot));
        while (!m_arpack_solver_nonsym->ArnoldiBasisFound()) {
            mpi::barrier();
            // Calling ARPACK FORTRAN code. Almost all work needed to
            // find an Arnoldi basis is performed by TakeStep.
            m_arpack_solver_nonsym->TakeStep();
            if ((m_arpack_solver_nonsym->GetIdo() == 1) || (m_arpack_solver_nonsym->GetIdo() == -1)) {
                mat.parallel_multiply(m_arpack_solver_nonsym->GetVector(),
                                      mat.m_nrow, m_arpack_solver_nonsym->PutVector());
            }
        }
        // Finding eigenvalues and eigenvectors.
        m_arpack_solver_nonsym->FindEigenvectors();
        return true;
    }

};


#endif //M7_ARNOLDISOLVER_H
