//
// Created by rja on 03/11/22.
//

#include "Rdms.h"
#include "SpinfreeRdm.h"
#include "Bilinears.h"

Rdms::exsig_to_rdms_t Rdms::make_exsig_to_rdms() const {
    exsig_to_rdms_t exsig_to_rdms;
    for (const auto& rdm: m_rdms) {
        auto contr_cast = dynamic_cast<const ContractedRdm*>(rdm.get());
        // pointer is only non-null (i.e. bool true) if the cast was successful
        auto max_contrib_exsig = contr_cast ? contr_cast->m_max_contrib_exsig : rdm->m_ranksig;
        auto nfrm_cre = max_contrib_exsig.nfrm_cre();
        auto nfrm_ann = max_contrib_exsig.nfrm_ann();
        while (nfrm_cre != ~0ul && nfrm_ann != ~0ul) {
            auto nbos_cre = max_contrib_exsig.nbos_cre();
            auto nbos_ann = max_contrib_exsig.nbos_ann();
            while (nbos_cre != ~0ul && nbos_ann != ~0ul) {
                OpSig exsig ({nfrm_cre, nfrm_ann}, {nbos_cre, nbos_ann});
                DEBUG_ASSERT_NE(exsig, opsig::c_invalid, "exsig OOB");
                exsig_to_rdms[exsig].push_front(rdm.get());
                --nbos_cre;
                --nbos_ann;
            }
            --nfrm_cre;
            --nfrm_ann;
        }
    }
    return exsig_to_rdms;
}

wf_comp_t Rdms::contrib_norm(uint_t iroot) const {
    const auto ipart1 = m_wf.m_format.flatten(iroot, 0);
    const auto ipart2 = m_wf.nreplica()==1 ? ipart1 : m_wf.m_format.flatten(iroot, 1);
    const auto& l2_norm_squares = m_wf.m_stats.m_l2_norm_square.total();
    return std::sqrt(l2_norm_squares[ipart1] * l2_norm_squares[ipart2]);
}

Rdms::Rdms(const conf::Rdms& opts, const wf::Vectors& wf, const Epoch& accum_epoch) :
        m_opts(opts), m_wf(wf), m_spinfree(opts.m_spinfree), m_work_conns(m_wf.m_sector.size()),
        m_work_com_ops(m_wf.m_sector.size()), m_accum_epoch(accum_epoch) {
    for (const auto& ranksig: bilinears::parse_exsigs(opts.m_ranks)) {
        REQUIRE_FALSE(m_pure_rdms[ranksig], "No RDM rank should appear more than once in the specification");
        REQUIRE_NE(ranksig, opsig::c_invalid, "invalid RDM rank signature (perhaps too many operators for OpSig object)");
        REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
        REQUIRE_TRUE(ranksig.conserves_nfrm(), "fermion non-conserving RDMs are not yet supported");
        REQUIRE_LE(ranksig.nbos(), 1ul, "RDMs with more than one boson operator are not yet supported");
        m_rdms.emplace_front(ptr::smart::make_poly_unique<Rdm, PureRdm>(opts, ranksig, m_wf.m_sector, 1ul));
        auto rdm_ptr = m_rdms.front().get();
        m_pure_rdms[ranksig] = rdm_ptr;
    }
    if (opts.m_fock_4rdm.m_enabled) {
        logging::info("Loading generalized Fock matrix for accumulation of its contraction with the 4RDM");
        FockMatrix fock(m_wf.m_sector.m_frm.m_basis.m_nsite, opts.m_fock_4rdm.m_fock_path);
        const auto screen_thresh = opts.m_fock_4rdm.m_screen_thresh.m_value;
        if (screen_thresh != 0.0) {
            logging::info("Screening Fock matrix element with magnitude lower than {:.2g} to zero", screen_thresh);
            fock.screen(screen_thresh);
        }
        const auto diag = fock.is_diagonal();
        logging::info("The given Fock matrix was found to be {}diagonal", (diag ? "" : "non-"));

        if (!diag) m_rdms.emplace_front(ptr::smart::make_poly_unique<Rdm, NonDiagFockRdm4>(opts, fock, m_wf.m_sector, 1ul));
        else m_rdms.emplace_front(ptr::smart::make_poly_unique<Rdm, DiagFockRdm4>(opts, fock, m_wf.m_sector, 1ul));
    }
    m_exsig_to_rdms = make_exsig_to_rdms();

    m_total_norm.m_local = 0.0;

    strv_t exsigs;
    for (uint_t iexsig=0ul; iexsig < opsig::c_ndistinct; ++iexsig){
        const OpSig exsig(iexsig);
        if (takes_contribs_from(exsig)) exsigs.push_back(exsig.to_string());
    }
    if (!exsigs.empty())
        logging::info("Excitation signatures contributing to the sampled RDMs: {}", convert::to_string(exsigs));

    if (m_opts.m_stoch_thresh_mag > 0.0) {
        const uint_t seed = 123;
        const uint_t block_size = 10000;
        m_stoch_thresh_prng = ptr::smart::make_unique<PRNG>(seed, block_size);
    }
    for (auto& rdm: m_rdms) {
        if (rdm->m_stoch_thresh_contribs) {
            logging::info("Stochastically applying standard-normalized threshold {:.2g} on contributions to RDM {}",
                          rdm->name(), m_opts.m_stoch_thresh_mag.m_value);
            REQUIRE_TRUE_ALL(m_stoch_thresh_prng.get(),
                             "Stochastic thresholding of RDM contribs requires non-zero specified positive magnitude");
        }
        if (rdm->m_neglect_tiny_contribs) {
            logging::info("Neglecting all |CiCj| contributions to RDM {} below magnitude {:.2g} in standard normalization",
                          rdm->name(), m_opts.m_stoch_thresh_mag.m_value);
        }
    }

}

Rdms::operator bool() const {
    return !m_rdms.empty();
}

bool Rdms::takes_contribs_from(OpSig exsig) const {
    return (exsig != opsig::c_invalid) && !m_exsig_to_rdms[exsig].empty();
}

void Rdms::make_contribs(const Mbf& src_onv, const conn::Mbf& conn, const com_ops::Mbf& com, const wf_t& contrib) {
    auto exsig = conn.exsig();
    if (exsig == opsig::c_invalid) return;
    if (!exsig) m_total_norm.m_local+=contrib;
    for (auto& rdm: m_exsig_to_rdms[exsig]) {
        auto rdm_contrib = contrib;
        if (rdm->m_neglect_tiny_contribs){
            const auto norm = this->contrib_norm(0);
            // skip this contribution if it is tiny
            if (std::abs(contrib / norm) < m_opts.m_neglect_contribs_mag) continue;
        }
        if (rdm->m_stoch_thresh_contribs){
            const auto norm = this->contrib_norm(0);
            const auto rdm_contrib = m_stoch_thresh_prng->stochastic_threshold(contrib / norm, m_opts.m_stoch_thresh_mag.m_value);
            // skip this contribution if it was stochastically rounded to zero
            if (rdm_contrib == 0.0) continue;
        }
        rdm->make_contribs(src_onv, conn, com, rdm_contrib);
    }
}

void Rdms::make_contribs(const Mbf& src_onv, const Mbf& dst_onv, const wf_t& contrib) {
    m_work_conns[src_onv].connect(src_onv, dst_onv, m_work_com_ops[src_onv]);
    make_contribs(src_onv, m_work_conns[src_onv], m_work_com_ops[src_onv], contrib);
}

bool Rdms::all_stores_empty() const {
    for (auto& rdm: m_rdms) if (!rdm->m_store.empty()) return false;
    return true;
}

void Rdms::end_cycle() {
    for (auto& rdm: m_rdms) rdm->end_cycle();
    m_total_norm.all_sum();
}

bool Rdms::is_energy_sufficient(const Hamiltonian& ham) const {
    if (ham.m_bos.m_contribs_0011.is_nonzero(opsig::c_zero)){
        if (!m_pure_rdms[opsig::c_0011]) return false;
    }
    if (ham.m_frmbos.m_contribs_1101.is_nonzero(opsig::c_1101)){
        if (!m_pure_rdms[opsig::c_1101]) return false;
    }
    if (ham.m_frmbos.m_contribs_1110.is_nonzero(opsig::c_1110)){
        if (!m_pure_rdms[opsig::c_1110]) return false;
    }
    if (!m_pure_rdms[opsig::c_doub]) return false;
    return true;
}


ham_comp_t Rdms::get_energy(const FrmHam& ham) const {
    if (!ham) return 0.0;
    auto& rdm2 = m_pure_rdms[opsig::c_doub];
    REQUIRE_TRUE_ALL(rdm2!=nullptr, "cannot compute energy without the 2RDM");
    ham_t e1 = 0.0;
    ham_t e2 = 0.0;
    wf_t trace = 0.0;
    auto& row = rdm2->m_store.m_row;

    for (row.restart(); row; ++row){
        const uint_t i=row.m_inds.m_frm.m_cre[0];
        const uint_t j=row.m_inds.m_frm.m_cre[1];
        DEBUG_ASSERT_LT(i, j, "spin orbital creation indices should be ordered");
        const uint_t k=row.m_inds.m_frm.m_ann[0];
        const uint_t l=row.m_inds.m_frm.m_ann[1];
        DEBUG_ASSERT_LT(k, l, "spin orbital annihilation indices should be ordered");
        const auto rdm_element = row.m_values[0];
        e2 += rdm_element*ham.get_coeff_2200(i, j, k, l);
        /*
         * signs of the following contributions come from the Fermi phase of bringing the like-valued creation and
         * annihilation operators together to act on the ket, leading to a factor of (nelec - 1), accounted for
         * after the loop.
         */
        if (i == k) e1 += rdm_element*ham.get_coeff_1100(j,l);
        if (j == l) e1 += rdm_element*ham.get_coeff_1100(i,k);
        if (i == l) e1 -= rdm_element*ham.get_coeff_1100(j,k);
        if (j == k) e1 -= rdm_element*ham.get_coeff_1100(i,l);
        if ((i==k) && (j==l)) trace+=rdm_element;
    }
    const uint_t nelec = rdm2->m_sector.m_frm.m_elecs;
    // scale the one-body contribution by the number of two-body contributions
    e1 /= nelec-1;
    e1 = mpi::all_sum(e1);
    e2 = mpi::all_sum(e2);
    trace = mpi::all_sum(trace);
    DEBUG_ASSERT_GT(std::abs(trace), 1e-14, "RDM trace should be non-zero");
    const auto norm = arith::real(trace) / integer::nspair(nelec);
    if (!rdm2->approx_contribs())
        REQUIRE_NEAR_EQ(norm, m_total_norm.m_reduced.real(),
                      "2RDM norm should match total of sampled diagonal contributions");
    REQUIRE_NEAR_ZERO(m_total_norm.m_reduced.imag(), "2RDM norm should be purely real");
    return arith::real(ham.m_e_core) + (arith::real(e1) + arith::real(e2)) / norm;
}

ham_comp_t Rdms::get_energy(const FrmBosHam& /*ham*/, OpSig /*exsig*/) const {
    return 0.0;
    // TODO: update for new HamOpTerm partitioning
#if 0
    if (!ham) return 0.0;
    auto& rdm = m_rdms[exsig];
    REQUIRE_TRUE_ALL(exsig::decode_nbos(exsig)==1,
                     "currently only supporting linear boson operators in the ladder term");
    REQUIRE_TRUE_ALL(rdm!=nullptr, "cannot compute energy without the "+exsig::to_string(exsig)+"-RDM");
    ham_t e_uncoupled = 0.0; // 0001 and 0010
    ham_t e_coupled = 0.0; // 1101 and 1110
    auto& row = rdm->m_store.m_row;
    bool cre = exsig::decode_nbos_cre(exsig);

    for (row.restart(); row.in_range(); row.step()){
        const uint_t p=row.m_inds.m_frm.m_cre[0];
        const uint_t q=row.m_inds.m_frm.m_ann[0];
        const uint_t n=cre ? row.m_inds.m_bos.m_cre[0] : row.m_inds.m_bos.m_ann[0];

        const auto rdm_element = row.m_values[0];
        /*
         * default definition is the creation ladder operator, so the fermion indices must be interchanged if computing
         * the energy of the boson-annihilating RDM
         */
        const auto coeff = cre ? ham.get_coeff_1110(n, p, q) : ham.get_coeff_1101(n, p, q);
        e_coupled += rdm_element * coeff;

        if (p == q) e_uncoupled += rdm_element*ham.get_coeff_0001(n);
    }
    e_uncoupled = mpi::all_sum(e_uncoupled);
    e_uncoupled/=nelec;
    e_coupled = mpi::all_sum(e_coupled);
    auto e = (e_uncoupled + e_coupled) / m_total_norm.m_reduced;
    REQUIRE_NEAR_EQ(dtype::imag(e), 0.0, 1e-12, "energy should be purely real");
    return dtype::real(e);
#endif
}

ham_comp_t Rdms::get_energy(const BosHam& ham) const {
    if (!ham) return 0.0;
    auto& rdm = m_pure_rdms[opsig::c_0011];
    REQUIRE_TRUE_ALL(rdm!=nullptr, "cannot compute energy without the 0011-RDM");
    ham_t e = 0.0;
    auto& row = rdm->m_store.m_row;

    for (row.restart(); row; ++row){
        const uint_t n=row.m_inds.m_bos.m_cre[0];
        const uint_t m=row.m_inds.m_bos.m_ann[0];
        REQUIRE_EQ(n, m, "0011-RDM should currently only take 0000-exsig contributions");
        const auto rdm_element = row.m_values[0];
        e += rdm_element*ham.get_coeff_0011(n, m);
    }
    e = mpi::all_sum(e) / m_total_norm.m_reduced;
    REQUIRE_TRUE(fptol::near_real(e), "energy should be purely real")
    return arith::real(e);
}

void Rdms::save(const hdf5::NodeWriter& parent) {
    if (!m_accum_epoch) {
        logging::warn("MAE accumulation epoch was not reached in this calculation: omitting RDM save");
        return;
    }
    {
        // unnormalized RDMs, suitable for restarts
        hdf5::GroupWriter gw(parent, "archive");
        hdf5::DatasetSaver::save_scalar(gw, "norm", wf_t(m_total_norm.m_reduced));
        for (const auto& rdm: m_rdms) rdm->save(gw);
    }

    if (m_spinfree) {
        /*
         * create and save the spinfree versions of all RDMs and intermediates
         */
        hdf5::GroupWriter gw(parent, "spinfree");
        for (const auto& rdm: m_rdms) SpinFreeRdm(*rdm, m_total_norm.m_reduced).save(gw);
    }
}

void Rdms::save() {
    if (!m_accum_epoch) {
        logging::warn("MAE accumulation epoch was not reached in this calculation: omitting RDM save");
        return;
    }
    REQUIRE_TRUE(m_opts.m_save.m_enabled, "save() called on Rdms object but saving was not enabled");
    save(hdf5::FileWriter(m_opts.m_save.m_path));
}

