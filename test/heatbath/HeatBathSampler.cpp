//
// Created by rja on 29/02/2020.
//

#include <gtest/gtest.h>
#include <src/sample/PRNG.h>
#include <src/defs.h>
#include <src/heatbath/HeatBathSampler.h>
#include <src/hamiltonian/AbInitioHamiltonian.h>
#include <src/heatbath/DeterminantSampler.h>

TEST(HeatBathSampler, AllExcitsGenerated) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    HeatBathSampler heat_bath_sampler(ham);
    auto source_det = ham.guess_reference(0);
    auto det_sampler = DeterminantSampler(heat_bath_sampler, source_det);
    PRNG prng(18);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0/nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    //all_connections.print();
    defs::inds frequencies(all_connections.high_water_mark(), 0ul);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        auto excit = det_sampler.draw(prng);
        if (!excit.m_single.is_null()) {
            auto det = excit.m_single.get_connection();
            //det.print();
            irow = all_connections.lookup(det);
            if (irow==~0ul) ASSERT_EQ(det.nexcit(source_det), 1);
            else frequencies[irow]++;
            frequencies[irow]++;
        }
        if (!excit.m_double.is_null()) {
            auto det = excit.m_double.get_connection();
            //det.print();
            irow = all_connections.lookup(det);
            if (irow==~0ul) ASSERT_EQ(det.nexcit(source_det), 2);
            else frequencies[irow]++;
        }
        if (std::all_of(frequencies.begin(), frequencies.end(), [](size_t i){return i>0;})){
            break;
        }
    }
    ASSERT_TRUE(std::all_of(frequencies.begin(), frequencies.end(), [](size_t i){return i>0;}));
}

TEST(HeatBathSampler, Unbiased) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    HeatBathSampler heat_bath_sampler(ham);
    auto source_det = ham.guess_reference(0);
    auto det_sampler = DeterminantSampler(heat_bath_sampler, source_det);

    PRNG prng(18);

    const size_t nattempt = 1e6;
    const defs::ham_comp_t eps = 1e-2;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    //all_connections.print();
    defs::inds frequencies(all_connections.high_water_mark(), 0ul);
    std::vector<defs::prob_t> probabilities(all_connections.high_water_mark(), 0);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        auto excit = det_sampler.draw(prng);
        if (!excit.m_single.is_null()) {
            auto det = excit.m_single.get_connection();
            //det.print();
            irow = all_connections.lookup(det);
            if (irow==~0ul) ASSERT_EQ(det.nexcit(source_det), 1);
            else {
                if (frequencies[irow]==0) probabilities[irow] = excit.m_single.m_prob;
                frequencies[irow]++;
            }
        }
        if (!excit.m_double.is_null()) {
            auto det = excit.m_double.get_connection();
            //det.print();
            irow = all_connections.lookup(det);
            if (irow==~0ul) ASSERT_EQ(det.nexcit(source_det), 2);
            else {
                if (frequencies[irow]==0) {
                    probabilities[irow] = excit.m_double.m_prob;
                    ASSERT_FALSE(consts::float_is_zero(probabilities[irow]));
                }
                frequencies[irow]++;
            }
        }
    }
    for (size_t i=0ul; i<frequencies.size(); ++i){
        std::cout << probabilities[i]/(frequencies.size()*frequencies[i]/(defs::prob_t)nattempt) <<std::endl;
    }
    //ASSERT_TRUE(std::all_of(frequencies.begin(), frequencies.end(), [](size_t i){return i>0;}));
}