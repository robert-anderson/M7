//
// Created by rja on 09/05/2020.
//

#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include "gtest/gtest.h"
#include "src/core/excitgen/HeatBathDoubles.h"
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"

namespace heat_bath_doubles_test {
    struct TestTableSpec : MappedTable<fields::FermionOnv> {
        size_t m_nattempt;
        fields::FermionOnv m_onv;
        fields::Numbers<size_t, 1> m_frequency;
        fields::Numbers<defs::prob_t, 1> m_weight;

        TestTableSpec(const views::FermionOnv &src_fonv, size_t nattempt) :
                MappedTable<fields::FermionOnv>(m_onv, 1000),
                m_nattempt(nattempt),
                m_onv(this, {src_fonv.nsite()}, "occupation number vector"),
                m_frequency(this, "number of times the ONV was drawn", nattempt),
                m_weight(this, "cumulative reciprocal probability", nattempt) {}
    };

    struct TestTable : BufferedTable<TestTableSpec> {
        TestTable(const FermionHamiltonian &ham, const views::FermionOnv &src_fonv, size_t nattempt) :
                BufferedTable<TestTableSpec>("Excitation generator testing mapped table", src_fonv, nattempt) {
            /*
             * table will expand dynamically
             */
            expand(10);
            HamiltonianDoubleConnectionEnumerator enumerator(ham, src_fonv);
            elements::FermionOnv dst_fonv(ham.nsite());
            MatrixElement<defs::ham_t> matel(src_fonv);
            while (enumerator.next(matel)) {
                matel.aconn.apply(src_fonv, dst_fonv);
                if (is_full()) expand_by_factor(1);
                insert(dst_fonv);
            }
        }

        bool all_have_at_least_one() const {
            for (size_t irow = 0; irow < m_hwm; ++irow)
                for (size_t iattempt = 0ul; iattempt < m_nattempt; ++iattempt)
                    if (!m_frequency(irow, iattempt)) return false;
            return true;
        }

        size_t errors_decreasing(size_t freq_cutoff, defs::inds ndraws) const {
            /*
             * returns row index of first counterexample
             */
            auto any_lt_cutoff = [&](size_t irow) {
                for (size_t iattempt = 0ul; iattempt < m_nattempt; ++iattempt)
                    if (m_frequency(irow, iattempt) < freq_cutoff) return true;
                return false;
            };

            for (size_t irow = 0; irow < m_hwm; ++irow) {
                // insufficient draws of this onv to be considered
                if (any_lt_cutoff(irow)) continue;
                for (size_t iattempt = 1ul; iattempt < m_nattempt; ++iattempt) {
                    if (
                            std::abs(1.0 - m_weight(irow, iattempt)/ndraws[iattempt]) >
                            std::abs(1.0 - m_weight(irow, iattempt - 1)/ndraws[iattempt]))
                        return irow;
                }
            }
            return ~0ul;
        }
    };
}


TEST(HeatBathDoubles, UnbiasedExcitsFromHFDeterminantRealSchroedinger) {
    Hamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false, 0, 0, 0);
    ASSERT_TRUE(ham.spin_conserving());
    PRNG prng(14, 1000000);
    HeatBathDoubles pchb(&ham, prng);

    elements::FermionOnv src_fonv(ham.nsite());
    elements::FermionOnv dst_fonv(ham.nsite());
    defs::inds occ_inds = {0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17};
    src_fonv.set(occ_inds);

    const size_t nattempt = 2;
    heat_bath_doubles_test::TestTable table(ham, src_fonv, nattempt);

    OccupiedOrbitals occ(src_fonv);
    VacantOrbitals vac(src_fonv);
    defs::prob_t prob;
    defs::ham_t helem;
    conn::AsFermionOnv aconn(src_fonv);

    const defs::inds ndraws = {1ul << 23, 1ul << 27};
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        const auto ndraw = ndraws[iattempt];
        std::cout << "ndraw: " << ndraw << std::endl;
        for (size_t idraw = 0ul; idraw < ndraw; ++idraw) {
            if (pchb.draw(src_fonv, dst_fonv, occ, vac, prob, helem, aconn)) {
                auto irow = *table[dst_fonv];
                table.m_frequency(irow, iattempt)++;
                table.m_weight(irow, iattempt) += (1.0 / prob);
            }
        }
    }

    ASSERT_TRUE(table.all_have_at_least_one());
    size_t irow = table.errors_decreasing(1e4, ndraws);
    if (irow != ~0ul) {
        std::cout << src_fonv.to_string() << std::endl;
        std::cout << table.m_onv(irow).to_string() << std::endl;
        std::cout << irow << " " <<
                  std::abs(1.0 - table.m_weight(irow, 0)) << " " <<
                  std::abs(1.0 - table.m_weight(irow, 1)) << std::endl;
        std::cout << irow << " " <<
                  table.m_frequency(irow, 0) << " " <<
                  table.m_frequency(irow, 1) << std::endl;
    }
    ASSERT_TRUE(table.errors_decreasing(1e4, ndraws) == ~0ul);
}

#if 0

bool excit_gen_tester(ExcitationGenerator &exgen, const FermionOnv &src_det, size_t ndraw, size_t pc_freq_thresh) {
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
                MappedList("test excitation connection list", m_determinant, nbucket),
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
    AntisymFermionOnvConnection anticonn(src_det);
    FermionOnv work_det(src_det);

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
            size_t irow = connection_list.lookup_irow(work_det);
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


TEST(HeatBathDoubles, UnbiasedExcitsFromHFDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    PRNG prng(18, 1e4);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    defs::inds occ_inds = {0, 1, 4, 5};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathDoubles, UnbiasedExcitsFromExcitedDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    PRNG prng(15, 10000);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    source_det.set(defs::inds{1, 4, 6, 7});
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathDoubles, UnbiasedExcitsFromSpinnedDeterminantComplex4c) {
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();

    AbInitioHamiltonian ham(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_FALSE(ham.spin_conserving());
    PRNG prng(15, 10000);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    source_det.set(defs::inds{1, 5, 6, 7});
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathDoubles, UnbiasedExcitsFromHFDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    PRNG prng(15, 10000);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    defs::inds occ_inds = {0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathDoubles, UnbiasedExcitsFromExcitedDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    PRNG prng(15, 10000);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    defs::inds occ_inds = {1, 4, 5, 7, 8, 10, 12, 14, 15, 16, 20, 21};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}

TEST(HeatBathDoubles, UnbiasedExcitsFromSpinnedDeterminantRealSchroedinger) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_Cr2_12o12e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    PRNG prng(15, 100000);
    HeatBathDoubles pchb(&ham, prng);

    FermionOnv source_det(ham.nsite());
    defs::inds occ_inds = {1, 4, 5, 7, 8, 9, 10, 14, 15, 16, 20, 21};
    source_det.set(occ_inds);
    ASSERT_TRUE(excit_gen_tester(pchb, source_det, 1e8, 1e6));
}
#endif
