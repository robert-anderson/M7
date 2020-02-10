//
// Created by Robert John Anderson on 2020-02-05.
//

#include <gtest/gtest.h>

#include "../../src/fermion/Determinant.h"
#include "../../src/enumerators/BitfieldEnumerator.h"
#include "../../src/parallel/MPIWrapper.h"


TEST(Determinant, Phase){

    Determinant bra(30);
    Determinant ket(30);
    bra.set(defs::inds{0, 1, 2, 3, 4, 5});
    ket.set(defs::inds{0, 1, 2, 3, 5, 8});

    std::cout << std::endl;
    bra.print();
    ket.print();
    std::cout << std::endl;
    std::cout << bra.phase(ket) << std::endl;
    std::cout << ket.phase(bra) << std::endl;
    MPIWrapper mpi_wrapper;
    ASSERT_LT(mpi_wrapper.irank(), mpi_wrapper.nrank());
}