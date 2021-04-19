//
// Created by rja on 12/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/excitgen/BosonExcitationGenerator.h"

#if 0
namespace boson_coupling_samplers_test {
    struct TestTable : MappedTable<fields::Onv<1>> {
        fields::Onv<1> m_onv;
        fields::Number<size_t> m_frequency;
        fields::Number<defs::prob_t> m_weight;

        TestTable(size_t nsite) :
                MappedTable<fields::Onv<1>>(m_onv, 1000),
                m_onv(this, nsite, "occupation number vector"),
                m_frequency(this, "number of times the ONV was drawn"),
                m_weight(this, "cumulative reciprocal probability") {}


        bool test_boson_gen(elements::FermiBosOnv &src_onv, size_t ndraw, size_t nboson_max) {
            size_t nsite = src_onv.m_frm.nsite();
            BufferedTable<boson_coupling_samplers_test::TestTable> bt("Excit gen tester", nsite);

            populate_table(bt, src_onv, ndraw, nboson_max);

            bool result = all_have_at_least_1(src_onv, nboson_max, bt.m_hwm);
            result &= errors_decreasing(src_onv, ndraw, nboson_max);
            result &= is_uniform(bt.m_frequency, ndraw);
            return result;
        }

        void populate_table(BufferedTable<boson_coupling_samplers_test::TestTable> &bt,
                            elements::FermiBosOnv &src_onv,
                            size_t ndraw,
                            size_t nboson_max) {
            size_t nsite = src_onv.m_frm.nsite();
            elements::FermiBosOnv dst_onv(nsite);
            PRNG prng = PRNG(18, 1e4);
            FermiBosHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false, nboson_max, 1.0, 0.5);
            BosonExcitationGenerator exgen(&ham, prng, nboson_max);
            OccupiedOrbitals occ_orbs(src_onv.m_frm);
            VacantOrbitals vac_orbs(src_onv.m_frm);
            conn::Antisym<1> aconn(src_onv);
            defs::prob_t prob;
            defs::ham_t helem;
            for (size_t idraw = 0ul; idraw < ndraw; ++idraw) {
                exgen.draw(src_onv, dst_onv, occ_orbs, vac_orbs, prob, helem, aconn);
                auto irow = *bt[dst_onv];
                if (irow == ~0ul) {
                    bt.expand(1);
                    irow = bt.insert(dst_onv);
                }
                bt.m_frequency(irow)++;
                bt.m_weight(irow) += 1.0 / prob;
            }
        }

        bool all_have_at_least_1(elements::FermiBosOnv &src_onv, size_t nboson_max, size_t num_generated) {
            // check correct no. of excitations is generated.
            size_t num_expected = 0;
            for (size_t imode = 0; imode < src_onv.m_bos.nmode(); ++imode) {
                if(src_onv.m_frm.get(0, imode) or src_onv.m_frm.get(1, imode)){
                    num_expected++;
                    num_expected += (src_onv.m_bos(imode) != 0 and src_onv.m_bos(imode) != nboson_max);
                }
            }
            return num_expected == num_generated;
        }

        bool
        errors_decreasing(elements::FermiBosOnv &src_onv, size_t ndraw_init, size_t nboson_max, size_t ntries = 3) {
            BufferedTable<boson_coupling_samplers_test::TestTable> bt("Error tester", src_onv.m_frm.nsite());
            std::vector<defs::prob_t> vars(2 * src_onv.m_bos.nmode(), 1000);
            for (size_t itry = 0; itry < ntries; ++itry) {
                size_t ndraw = (size_t) ndraw_init * std::pow(100, itry);
                populate_table(bt, src_onv, ndraw, nboson_max);
                for (size_t irow = 0; irow < bt.m_hwm; ++irow) {
                    auto this_var = std::abs(1.0 - bt.m_weight(irow) / ndraw);
                    if (this_var > vars[irow]) return false;
                    vars[irow] = this_var;
                }
            }
            return true;
        }

        bool is_uniform(fields::Number<size_t> &freqs, size_t ndraw, double tol = 2.0) {
            size_t min_freq = ~0ul;
            size_t max_freq = 0;
            for (size_t ifreq = 0; ifreq < freqs.m_field.m_nelement; ++ifreq) {
                auto this_freq = freqs(ifreq);
                min_freq = this_freq < min_freq ? this_freq : min_freq;
                max_freq = this_freq > max_freq ? this_freq : max_freq;
            }
            return std::abs(max_freq - (double) min_freq) < std::sqrt(ndraw / 12.0) / tol;
        }
    };
}



TEST(BosonExcitationGenerator, SingleOnvTest){
    const size_t nsite = 6;
    const size_t nboson_max = 3;

    const size_t ndraw = 100000;
    elements::FermiBosOnv src_onv(nsite);
    src_onv = {{0, 4, 6, 11}, {1, 0, 0, 1, 3, 2}};

    boson_coupling_samplers_test::TestTable tt(nsite);

    tt.test_boson_gen(src_onv, ndraw, nboson_max);

}
#endif //ENABLE_BOSONS