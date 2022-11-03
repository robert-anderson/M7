//
// Created by rja on 03/11/22.
//

#include "Rdms.h"
#include "SpinfreeRdm.h"

std::array<uintv_t, exsig::c_ndistinct> Rdms::make_exsig_ranks() const {
    std::array<uintv_t, exsig::c_ndistinct> exsig_ranks;
    for (const auto& ranksig: m_rdm_ranksigs) {
        auto nfrm_cre = decode_nfrm_cre(ranksig);
        auto nfrm_ann = decode_nfrm_cre(ranksig);
        while (nfrm_cre != ~0ul && nfrm_ann != ~0ul) {
            auto nbos_cre = decode_nbos_cre(ranksig);
            auto nbos_ann = decode_nbos_ann(ranksig);
            while (nbos_cre != ~0ul && nbos_ann != ~0ul) {
                auto exsig = encode(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann);
                DEBUG_ASSERT_LT(exsig, exsig::c_ndistinct, "exsig OOB");
                exsig_ranks[exsig].push_back(ranksig);
                --nbos_cre;
                --nbos_ann;
            }
            --nfrm_cre;
            --nfrm_ann;
        }
    }
    return exsig_ranks;
}

Rdms::Rdms(const conf::Rdms& opts, uintv_t ranksigs, sys::Sector sector, const Epoch& accum_epoch) :
        Archivable("rdms", opts.m_archivable),
        m_spinfree(opts.m_spinfree), m_rdm_ranksigs(ranksigs), m_exsig_ranks(make_exsig_ranks()),
        m_work_conns(sector.size()), m_work_com_ops(sector.size()),
        m_explicit_ref_conns(opts.m_explicit_ref_conns), m_accum_epoch(accum_epoch) {
    for (const auto& ranksig: ranksigs) {
        REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
        REQUIRE_TRUE(conserves_nfrm(ranksig), "fermion non-conserving RDMs are not yet supported");
        REQUIRE_LE(decode_nbos_cre(ranksig), 1ul,
                   "RDMs with more than one boson creation operator are not yet supported");
        REQUIRE_LE(decode_nbos_ann(ranksig), 1ul,
                   "RDMs with more than one boson annihilation operator are not yet supported");
        REQUIRE_TRUE(m_rdms[ranksig] == nullptr, "No RDM rank should appear more than once in the specification");
        m_rdms[ranksig] = ptr::smart::make_unique<Rdm>(opts, ranksig, ranksig, sector, 1ul);
    }
    if (opts.m_fock_4rdm.m_enabled) {
        FockMatrix fock(sector.m_frm.size(), opts.m_fock_4rdm.m_fock_path);
        m_fock_rdm4 = ptr::smart::make_unique<FockRdm4>(opts, fock, sector, 1ul);
    }

    m_total_norm.m_local = 0.0;

    strv_t exsigs;
    for (uint_t exsig=0ul; exsig < exsig::c_ndistinct; ++exsig){
        if (takes_contribs_from(exsig)) exsigs.push_back(exsig::to_string(exsig));
    }
    logging::info("Excitation signatures contributing to the sampled RDMs: {}", convert::to_string(exsigs));

}

Rdms::operator bool() const {
    return !m_rdm_ranksigs.empty() || m_fock_rdm4;
}

bool Rdms::takes_contribs_from(uint_t exsig) const {
    if (exsig > exsig::c_ndistinct) return false;
//    if ((exsig==ex_quadruple) && m_fock_rdm4 && !m_fock_rdm4->m_diagonal) return true;
    if ((exsig==ex_quadruple) && m_fock_rdm4) return true;
    return !m_exsig_ranks[exsig].empty();
}

void Rdms::make_contribs(const Mbf& src_onv, const conn::Mbf& conn, const com_ops::Mbf& com, const wf_t& contrib) {
    auto exsig = conn.exsig();
    if (!exsig) m_total_norm.m_local+=contrib;
    for (auto ranksig: m_exsig_ranks[exsig]) {
        if (m_rdms[ranksig]) m_rdms[ranksig]->make_contribs(src_onv, conn, com, contrib);
    }
    if (m_fock_rdm4) {
        if (is_pure_frm(exsig) && decode_nfrm_cre(exsig) <= 4ul && decode_nfrm_ann(exsig) <= 4ul)
            m_fock_rdm4->make_contribs(src_onv, conn, com, contrib);
    }
}

void Rdms::make_contribs(const Mbf& src_onv, const Mbf& dst_onv, const wf_t& contrib) {
    m_work_conns[src_onv].connect(src_onv, dst_onv, m_work_com_ops[src_onv]);
    make_contribs(src_onv, m_work_conns[src_onv], m_work_com_ops[src_onv], contrib);
}

bool Rdms::all_stores_empty() const {
    for (auto& ranksig: m_rdm_ranksigs)
        if (!m_rdms[ranksig]->m_store.empty())
            return false;
    if (m_fock_rdm4) return m_fock_rdm4->m_store.empty();
    return true;
}

void Rdms::end_cycle() {
    for (auto& ranksig: m_rdm_ranksigs) m_rdms[ranksig]->end_cycle();
    if (m_fock_rdm4) m_fock_rdm4->end_cycle();
    m_total_norm.all_sum();
}

bool Rdms::is_energy_sufficient(const Hamiltonian& ham) const {
    if (ham.m_bos.m_contribs_0011.is_nonzero(0ul)){
        if (!m_rdms[exsig::ex_0011]) return false;
    }
    if (ham.m_frmbos.m_contribs_1101.is_nonzero(exsig::ex_1101)){
        if (!m_rdms[exsig::ex_1101]) return false;
    }
    if (ham.m_frmbos.m_contribs_1110.is_nonzero(exsig::ex_1110)){
        if (!m_rdms[exsig::ex_1110]) return false;
    }
    if (!m_rdms[exsig::ex_double]) return false;
    return true;
}


ham_comp_t Rdms::get_energy(const FrmHam& ham) const {
    if (!ham) return 0.0;
    auto& rdm2 = m_rdms[ex_2200];
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
    REQUIRE_NEARLY_EQ(norm, m_total_norm.m_reduced, "2RDM norm should match total of sampled diagonal contributions");
    return arith::real(ham.m_e_core) + (arith::real(e1) + arith::real(e2)) / norm;
}

ham_comp_t Rdms::get_energy(const FrmBosHam& /*ham*/, uint_t /*exsig*/) const {
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
    REQUIRE_NEARLY_EQ(dtype::imag(e), 0.0, 1e-12, "energy should be purely real");
    return dtype::real(e);
#endif
}

ham_comp_t Rdms::get_energy(const BosHam& ham) const {
    if (!ham) return 0.0;
    auto& rdm = m_rdms[exsig::ex_0011];
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

void Rdms::save_fn(const hdf5::NodeWriter& parent) {
    if (!m_accum_epoch) {
        logging::warn("MAE accumulation epoch was not reached in this calculation: omitting RDM save");
        return;
    }
    hdf5::GroupWriter gw(parent, "rdms");
    gw.write_data("norm", m_total_norm.m_reduced);
    for (const auto& i: m_rdm_ranksigs) {
        DEBUG_ASSERT_TRUE(m_rdms[i].get(), "active ranksig was not allocated!");
        m_rdms[i]->save(gw);
    }
    if (m_fock_rdm4) m_fock_rdm4->save(gw);

    if (m_spinfree) {
        /*
         * create and save the spinfree versions of all RDMs and intermediates
         */
        for (const auto& i: m_rdm_ranksigs) {
            DEBUG_ASSERT_TRUE(m_rdms[i].get(), "active ranksig was not allocated!");
            SpinFreeRdm(*m_rdms[i], m_total_norm.m_reduced).save(gw);
        }
        if (m_fock_rdm4) SpinFreeRdm(*m_fock_rdm4, m_total_norm.m_reduced).save(gw);
    }
}