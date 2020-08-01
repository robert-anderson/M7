//
// Created by Robert John Anderson on 2020-01-24.
//

#include <gtest/gtest.h>
#include <src/core/linalg/EigenSolver.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/linalg/DenseHamiltonian.h"

TEST(DenseHamiltonian, FciEnergyCheck4c) {
    if (!consts::is_complex<defs::ham_comp_t>()) GTEST_SKIP();
    DenseHamiltonian ham(AbInitioHamiltonian(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false));
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -14.40597603432, 1e-10));
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[1], -14.28883698406, 1e-10));
}

TEST(DenseHamiltonian, FciEnergyCheckRhf) {
    AbInitioHamiltonian ham_src(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    // compare the ground and first excited states to BAGEL's values
    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], -108.81138657563143, 1e-10));
}

TEST(DenseHamiltonian, h2o) {

    FcidumpFileReader<double> file_reader("/Users/robertjohnanderson/tmp/FCIDUMP.out", false);
    ASSERT_TRUE(file_reader.spin_resolved());
    AbInitioHamiltonian ham_src(file_reader);
    auto hf_energy = 0;//-98.73369952192;
    auto ci_energy = 0;//-98.81346551;
    auto ref = ham_src.guess_reference(0);
    //ref.zero();
    //ref.set("1001|1001");
    std::cout << ham_src.get_energy(ref)-hf_energy << std::endl;
    DenseHamiltonian ham(ham_src);
    auto solver = ham.diagonalize();
    std::cout << solver.m_evals[0]-ci_energy << std::endl;

    /*
    -98.7337
    -98.8571
     */


    //-8005.37440388

//    auto hf_energy = -76.07451146988;
//    auto ci_energy = -76.07535316;
//    auto bagel_ci_energy = -76.07535816;
//    std::cout << ci_energy-bagel_ci_energy << std::endl;
//
//    auto bagel_hf_energy_2zfit = -76.07449540;
//    auto bagel_hf_energy_3zfit = -76.07451639;
//    auto bagel_hf_energy_4zfit = -76.07452174;
//    std::cout << hf_energy-bagel_hf_energy_2zfit << std::endl;
//    std::cout << hf_energy-bagel_hf_energy_3zfit << std::endl;
//    std::cout << hf_energy-bagel_hf_energy_4zfit << std::endl;
//    FcidumpFileReader<double> file_reader("/Users/robertjohnanderson/tmp/FCIDUMP", true);
//    ASSERT_TRUE(file_reader.spin_resolved());
//    AbInitioHamiltonian ham_src(file_reader);
//    std::cout << ham_src.guess_reference(0).to_string() << std::endl;
//    std::cout << ham_src.get_energy(ham_src.guess_reference(0))-hf_energy;
//    DenseHamiltonian ham(ham_src);
//    auto solver = ham.diagonalize();
//    // compare the ground and first excited states to DIRAC's values
//    ASSERT_TRUE(consts::floats_nearly_equal(solver.m_evals[0], ci_energy, 1e-7));
}