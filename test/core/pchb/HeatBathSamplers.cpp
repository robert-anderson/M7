//
// Created by rja on 09/05/2020.
//

#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include "gtest/gtest.h"
#include "src/core/pchb/HeatBathSamplers.h"

bool excit_gen_tester(ExcitationGenerator &exgen, const Determinant &src_det, size_t ndraw, size_t pc_freq_thresh) {
    /*
     * given:
     *      ExcitationGenerator subclass instance
     *      source determinant
     *      number of draws "ndraw"
     *      1% threshold frequency "pc_freq_thresh"
     *
     * this function will perform ndraw singles and ndraw doubles draws from src_det,
     * storing the generation frequencies f1, and the cumulative weighted probabilities w1
     *
     * then a further n draws will be performed, accumulating a copy f2 of the frequencies f1
     * and a copy w2 of the weighted probabilities w1.
     *
     * for the test to pass, the dst_dets indexed i with a f2[i]>=pc_freq_thresh, must be have a
     * normalized w2[i] correct within 1%, and the normalized w2[i] must be closer to 1 than w1[i]
     *
     * for now, this is a fairly rough check, with weak criteria, and a more statistically rigorous
     * scheme could certainly be devised
     *
     * a successful test is indicated by a non-zero return value
     */

    defs::prob_t dn = ndraw;

    const defs::ham_comp_t eps = 100.0 / ndraw;

    struct ExcitConnectionList : public MappedList<DeterminantElement> {
        DeterminantField m_determinant;
        NumericField<defs::ham_t> m_helement;
        NumericField<size_t> m_frequency;
        NumericField<defs::prob_t> m_weight;

        ExcitConnectionList(size_t nsite, size_t nbucket) :
                MappedList(m_determinant, nbucket),
                m_determinant(this, 1, nsite), m_helement(this),
                m_frequency(this, 2), m_weight(this, 2) {}
    };
    size_t nsite = src_det.nsite();
    size_t nconn = integer_utils::combinatorial(2 * nsite, 4); // comfortable upper bound
    ExcitConnectionList connection_list(nsite, nconn);
    HamiltonianConnectionEnumerator enumerator(*exgen.ham(), src_det, eps);
    MatrixElement<defs::ham_t> matel(src_det);
    auto dst_det = src_det;
    while (enumerator.next(matel)) {
        matel.aconn.apply(src_det, dst_det);
        auto irow = connection_list.expand_push(dst_det);
        connection_list.m_helement(irow) = matel.element;
    }
    nconn = connection_list.high_water_mark(0);

    size_t nnull = 0ul;

    OccupiedOrbitals occ(src_det);
    VacantOrbitals vac(src_det);
    AntisymConnection anticonn(src_det);
    Determinant work_det(src_det);

    auto loop = [&](size_t ielement) {
        defs::prob_t prob;
        defs::ham_t helem;
        bool valid;
        for (size_t i = 0ul; i < 2 * ndraw; ++i) {
            if (i < ndraw) {
                valid = exgen.draw_double(src_det, work_det, occ, prob, helem, anticonn);
                if (valid) ASSERT(anticonn.nexcit() == 2)
            } else {
                valid = exgen.draw_single(src_det, work_det, occ, vac, prob, helem,anticonn);
                if (valid) ASSERT(anticonn.nexcit() == 1)
            }

            if (!valid) {
                ++nnull;
                continue;
            }
            ASSERT(src_det.nsetbit() == work_det.nsetbit())
            size_t irow = connection_list.lookup(work_det);
            if (irow != ~0ul) {
                connection_list.m_weight(irow, 0, ielement) += 1.0 / prob;
                connection_list.m_frequency(irow, 0, ielement) += 1;
            }
        }
        for (size_t irow = 0ul; irow < nconn; ++irow) std::cout << *connection_list.m_weight(irow, 0, ielement) << " ";
        std::cout << std::endl;
        for (size_t irow = 0ul; irow < nconn; ++irow)
            std::cout << *connection_list.m_frequency(irow, 0, ielement) << " ";
        std::cout << std::endl;
    };
    loop(0);
    loop(1);
    defs::prob_t tot_err1 = 0;
    defs::prob_t tot_err2 = 0;
    for (size_t i = 0ul; i < nconn; ++i) {
        auto f1 = *connection_list.m_frequency(i, 0, 0);
        auto f2 = *connection_list.m_frequency(i, 0, 1);
        f2 += f1;
        auto w1 = *connection_list.m_weight(i, 0, 0);
        auto w2 = *connection_list.m_weight(i, 0, 1);
        w2 += w1;
        // check that all accessible dst_dets have been generated at least once
        if (f2 == 0) return false;
        auto err1 = std::abs(w1 / (dn) - 1.0);
        auto err2 = std::abs(w2 / (2 * dn) - 1.0);
        tot_err1 += err1;
        tot_err2 += err2;
        if (f2 > pc_freq_thresh) {
            if (err2 > 0.01) {
                std::cout << err2 << " " << f2 << std::endl;
                return false;
            }
        }
    }
    return tot_err2 < tot_err1;
}


TEST(HeatBathSamplers, UnbiasedExcitsFromHFDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(18, 1e4);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    defs::inds occ_inds = {0, 1, 4, 5};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromExcitedDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();

    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(15, 10000);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    source_det.set(defs::inds{1, 4, 6, 7});
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromSpinnedDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();

    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    PRNG prng(15, 10000);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    source_det.set(defs::inds{1, 5, 6, 7});
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromHFDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP");
    PRNG prng(15, 10000);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    defs::inds occ_inds = {0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromExcitedDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP");
    PRNG prng(15, 10000);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    defs::inds occ_inds = {1, 4, 5, 7, 8, 10, 12, 14, 15, 16, 20, 21};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathSamplers, UnbiasedExcitsFromSpinnedDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP");
    PRNG prng(15, 100000);
    HeatBathSamplers pchb(&ham, prng);

    Determinant source_det(ham.nsite());
    defs::inds occ_inds = {1, 4, 5, 7, 8, 9, 10, 14, 15, 16, 20, 21};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}