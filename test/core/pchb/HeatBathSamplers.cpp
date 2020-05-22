//
// Created by rja on 09/05/2020.
//

#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include "gtest/gtest.h"
#include "src/core/pchb/HeatBathSamplers.h"

TEST(HeatBathSamplers, AllExcitsGeneratedFromHartreeFockDeterminantComplex4c) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PrivateStore<PRNG> prng(1, PRNG(18, 1e4));
    HeatBathSamplers pchb(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    OccupiedOrbitals occ(source_det);
    VacantOrbitals vac(source_det);
    AntisymConnection anticonn(source_det);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);

    size_t irow;
    defs::prob_t prob;
    defs::ham_t helem;
    size_t nnull=0ul;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        bool valid = pchb.draw_double(source_det, work_det, occ, prob, helem, anticonn);
        if (!valid){
            ++nnull;
            continue;
        }
        irow = all_connections.lookup(work_det);
        if (irow != ~0ul) frequencies[irow]++;
        if (std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; })) break;
    }
    ASSERT_EQ(nnull, 5239381);
    ASSERT_EQ(std::accumulate(frequencies.begin(), frequencies.end(), 0ul), 672869);
    ASSERT_TRUE(std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; }));
}

TEST(HeatBathSamplers, UnbiasedDoublesFromHartreeFockDeterminantComplex4c) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PrivateStore<PRNG> prng(1, PRNG(16, 3e4));
    HeatBathSamplers pchb(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    OccupiedOrbitals occ(source_det);
    VacantOrbitals vac(source_det);
    AntisymConnection anticonn(source_det);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);
    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    defs::prob_t prob;
    defs::ham_t helem;
    size_t nnull=0ul;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        bool valid = pchb.draw_double(source_det, work_det, occ, prob, helem, anticonn);
        if (!valid){
            ++nnull;
            continue;
        }
        ASSERT_EQ(source_det.nsetbit(), work_det.nsetbit());
        irow = all_connections.lookup(work_det);
        if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / (prob*nattempt);
        if (irow != ~0ul) frequencies[irow]++;
    }
    //utils::print(weighted_frequencies);
    //utils::print(frequencies);
    ASSERT_EQ(nnull, 88623452);
    ASSERT_TRUE(std::all_of(weighted_frequencies.cbegin(), weighted_frequencies.cend(), [](const defs::prob_t i) { return i > 0; }));
}


TEST(HeatBathSamplers, UnbiasedExcitsFromExcitedDeterminantComplex4c) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PrivateStore<PRNG> prng(1, PRNG(16, 3e4));
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    /*
     * arbitrary choice of source determinant
     */
    source_det.set(defs::inds{1, 4, 6, 7});
    Determinant work_det(ham.nsite());
    OccupiedOrbitals occ(source_det);
    VacantOrbitals vac(source_det);
    AntisymConnection anticonn(source_det);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);
    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    defs::prob_t prob;
    defs::ham_t helem;
    size_t nnull=0ul;
    bool valid;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        if (iattempt<nattempt/2)
            valid = pchb.draw_double(source_det, work_det, occ, prob, helem, anticonn);
        else
            valid = pchb.draw_single(source_det, work_det, occ, vac, prob, helem, anticonn);
        if (!valid){
            ++nnull;
            continue;
        }
        ASSERT_EQ(source_det.nsetbit(), work_det.nsetbit());
        irow = all_connections.lookup(work_det);
        if (irow != ~0ul) weighted_frequencies[irow] += 2.0 / (prob*nattempt);
        if (irow != ~0ul) frequencies[irow]++;
    }
    utils::print(weighted_frequencies);
    utils::print(frequencies);
    std::cout << nnull <<std::endl;
    ASSERT_EQ(nnull, 45562755);
    ASSERT_TRUE(std::all_of(weighted_frequencies.cbegin(), weighted_frequencies.cend(), [](const defs::prob_t i) { return i > 0; }));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromExcitedDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP");
    PrivateStore<PRNG> prng(1, PRNG(16, 3e4));
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    /*
     * arbitrary choice of source determinant
     */
    source_det.set(defs::inds{1, 4, 5, 7, 8, 10,   12, 14, 15, 16, 20, 21});
    Determinant work_det(ham.nsite());
    OccupiedOrbitals occ(source_det);
    VacantOrbitals vac(source_det);
    AntisymConnection anticonn(source_det);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);
    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    defs::prob_t prob;
    defs::ham_t helem;
    size_t nnull=0ul;
    bool valid;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        if (iattempt<nattempt/2)
            valid = pchb.draw_double(source_det, work_det, occ, prob, helem, anticonn);
        else
            valid = pchb.draw_single(source_det, work_det, occ, vac, prob, helem, anticonn);
        if (!valid){
            ++nnull;
            continue;
        }
        ASSERT_EQ(source_det.nsetbit(), work_det.nsetbit());
        irow = all_connections.lookup(work_det);
        if (irow != ~0ul) weighted_frequencies[irow] += 2.0 / (prob*nattempt);
        if (irow != ~0ul) frequencies[irow]++;
    }
    utils::print(weighted_frequencies);
    utils::print(frequencies);
    std::cout << nnull <<std::endl;
    ASSERT_EQ(nnull, 48899841);
    ASSERT_TRUE(std::all_of(weighted_frequencies.cbegin(), weighted_frequencies.cend(), [](const defs::prob_t i) { return i > 0; }));
}

