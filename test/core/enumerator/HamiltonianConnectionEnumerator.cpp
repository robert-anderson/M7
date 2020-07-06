//
// Created by rja on 06/07/2020.
//

#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include "gtest/gtest.h"
#include "src/core/enumerator/HamiltonianConnectionEnumerator.h"

TEST(HamiltonianConnectionEnumerator, Test){
    AbInitioHamiltonian ham(defs::assets_root+"/RHF_N2_6o6e/FCIDUMP");


    auto det = ham.guess_reference(0);
    det.zero();
    det.set({0, 1, 2, 6, 7, 8});
    HamiltonianConnectionEnumerator enumerator(ham, det);

    auto excited = det;
    MatrixElement<defs::ham_t> matrix_element(det);

    size_t count = 0;
    while (enumerator.next(matrix_element)){
        matrix_element.aconn.apply(det, excited);
        ASSERT_EQ(excited.nsetbit(), det.nsetbit());
        count++;
    }
    ASSERT_EQ(count, 47);
}
