//
// Created by rja on 29/02/2020.
//

#include <gtest/gtest.h>
#include <src/defs.h>
#include <src/core/heatbath/HeatBathSampler.h>
#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include <src/core/heatbath/DeterminantSampler.h>

TEST(HeatBathSampler, AllExcitsGeneratedFromHartreeFockDeterminantComplex4c) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e8;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw();

        // Brillouin's theorem:
        ASSERT_FALSE(det_sampler.single_generated());

        if (det_sampler.double_generated()) {
            det_sampler.get_double().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) frequencies[irow]++;
        }
        if (std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; })) break;
    }
    ASSERT_EQ(std::accumulate(frequencies.begin(), frequencies.end(), 0ul), 786387);
    ASSERT_TRUE(std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; }));
}

TEST(HeatBathSampler, AllExcitsGeneratedFromExcitedDeterminantComplex4c) {
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    Determinant source_det(ham.nsite());
    /*
     * arbitrary choice of source determinant
     */
    source_det.set(defs::inds{1, 4, 6, 7});
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e6;
    const defs::ham_comp_t eps = 100.0 / nattempt;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    defs::inds frequencies(all_connections.high_water_mark(0), 0ul);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw();
        if (det_sampler.single_generated()) {
            det_sampler.get_single().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) frequencies[irow]++;
        }
        if (det_sampler.double_generated()) {
            det_sampler.get_double().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) frequencies[irow]++;
        }
        if (std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; })) {
            break;
        }
    }
    ASSERT_EQ(std::accumulate(frequencies.begin(), frequencies.end(), 0ul), 2417);
    ASSERT_TRUE(std::all_of(frequencies.begin(), frequencies.end(), [](size_t i) { return i > 0; }));
}

TEST(HeatBathSampler, UnbiasedElecPairComplex4c) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e7;
    NdArray<defs::prob_t, 2> weighted_frequencies(ham.nsite() * 2, ham.nsite() * 2);

    size_t iattempt;
    defs::prob_t tol = 2e-3;
    auto predicate = [&iattempt, tol](defs::prob_t weight) {
        if (iattempt < 100) return false;
        return consts::float_is_zero(weight) ||
               consts::floats_nearly_equal(weight / (iattempt + 1), 1.0, tol);
    };

    size_t p, q;
    defs::prob_t prob;
    for (iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw_pq(p, q, prob);
        ASSERT_GT(prob, 0.0);
        ASSERT_LE(prob, 1.0);
        *weighted_frequencies.view(p, q) += 1.0 / prob;

        auto begin = weighted_frequencies.view(0, 0);
        if (std::all_of(begin, begin + weighted_frequencies.nelement(), predicate)) break;
    }

    auto begin = weighted_frequencies.view(0, 0);
    ASSERT_TRUE(std::all_of(begin, begin + weighted_frequencies.nelement(), predicate));
    ASSERT_EQ(iattempt, 5824562);
}

TEST(HeatBathSampler, UnbiasedElecPairAndFirstVirtualComplex4c) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e6;
    const auto nspinorb = ham.nsite() * 2;
    NdArray<defs::prob_t, 3> weighted_frequencies(nspinorb, nspinorb, nspinorb);

    size_t iattempt;
    defs::prob_t tol = 5e-3;
    auto predicate = [&iattempt, tol](defs::prob_t weight) {
        if (iattempt < 100) return false;
        return consts::float_is_zero(weight) ||
               consts::floats_nearly_equal(weight / (iattempt + 1), 1.0, tol);
    };

    size_t p, q, r;
    defs::prob_t prob;
    for (iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw_pqr(p, q, r, prob);
        if (r == ~0ul) continue;
        ASSERT_GT(prob, 0.0);
        ASSERT_LE(prob, 1.0);
        *weighted_frequencies.view(p, q, r) += 1.0 / prob;

        auto begin = weighted_frequencies.view(0, 0, 0);
        if (std::all_of(begin, begin + weighted_frequencies.nelement(), predicate)) break;
    }

    for (size_t p = 0ul; p < nspinorb; ++p) {
        for (size_t q = 0ul; q < nspinorb; ++q) {
            for (size_t r = 0ul; r < nspinorb; ++r) {
                if (*weighted_frequencies.view(p, q, r) > 0)
                    std::cout << p << " " << q << " " << r << " " <<
                              *weighted_frequencies.view(p, q, r) / (iattempt + 1) - 1.0
                              << std::endl;
            }
        }
    }
    /*
    auto begin = weighted_frequencies.view(0, 0, 0);
    ASSERT_TRUE(std::all_of(begin, begin + weighted_frequencies.nelement(), predicate));
    ASSERT_EQ(iattempt, 5824562);
    */

}

TEST(HeatBathSampler, UnbiasedFromHartreeFockDeterminantComplex4c) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */

    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e6;
    const defs::ham_comp_t eps = 1e-2;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw();
        if (det_sampler.single_generated()) {
            det_sampler.get_single().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_single_prob();
        }
        if (det_sampler.double_generated()) {
            det_sampler.get_double().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_double_prob();
        }
    }
    ASSERT_TRUE(std::all_of(weighted_frequencies.begin(), weighted_frequencies.end(), [](size_t i) { return i > 0; }));
}

TEST(HeatBathSampler, UnbiasedFromExcitedDeterminantComplex4c) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */

    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    Determinant source_det(ham.nsite());
    /*
     * arbitrary choice of source determinant
     */
    source_det.set(defs::inds{1, 4, 6, 7});
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e7;
    const defs::ham_comp_t eps = 1e-2;
    auto all_connections = ham.all_connections_of_det(source_det, eps);
    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw();
        if (det_sampler.single_generated()){
            det_sampler.get_single().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_single_prob();
        }
        if (det_sampler.double_generated()){
            det_sampler.get_double().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_double_prob();
        }
    }
    ASSERT_TRUE(std::all_of(weighted_frequencies.begin(), weighted_frequencies.end(), [](size_t i) { return i > 0; }));
}

TEST(HeatBathSampler, UnbiasedFromHartreeFockDeterminantRealSchroedinger) {
    /*
     * ensure that the ratio of generation frequency to proposal probability is uniform
     */
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSampler heat_bath_sampler(&ham, prng);
    auto source_det = ham.guess_reference(0);
    Determinant work_det(ham.nsite());
    auto &det_sampler = heat_bath_sampler.det_sampler->get(0);
    det_sampler.update(source_det);

    const size_t nattempt = 1e6;
    const defs::ham_comp_t eps = 1e-2;
    auto all_connections = ham.all_connections_of_det(source_det, eps);

    for (size_t irow = 0ul; irow<all_connections.nrow_per_segment(); ++irow)
        ASSERT_EQ(all_connections.determinant(irow).spin(), 0);

    std::vector<defs::prob_t> weighted_frequencies(all_connections.high_water_mark(0), 0);

    size_t irow;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        det_sampler.draw();
        if (det_sampler.single_generated()){
            det_sampler.get_single().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_single_prob();
        }
        if (det_sampler.double_generated()){
            det_sampler.get_double().apply(source_det, work_det);
            irow = all_connections.lookup(work_det);
            if (irow != ~0ul) weighted_frequencies[irow] += 1.0 / det_sampler.get_double_prob();
        }
    }
    ASSERT_TRUE(std::all_of(weighted_frequencies.begin(), weighted_frequencies.end(), [](size_t i) { return i > 0; }));
}