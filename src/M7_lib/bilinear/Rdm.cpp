//
// Created by rja on 11/08/2021.
//

#include "Rdm.h"

size_t Rdm::nrow_estimate(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann, const BasisData& bd) {
    double nrow = 1.0;
    nrow *= integer_utils::combinatorial(bd.m_frm.m_nspinorb, nfrm_cre);
    nrow *= integer_utils::combinatorial(bd.m_frm.m_nspinorb, nfrm_ann);
    nrow *= integer_utils::combinatorial(bd.m_bos.m_nmode, nbos_cre);
    nrow *= integer_utils::combinatorial(bd.m_bos.m_nmode, nbos_ann);
    nrow /= integer_utils::factorial(nfrm_cre + nfrm_ann);
    nrow /= integer_utils::factorial(nbos_cre + nbos_ann);
    return nrow;
}

size_t Rdm::nrow_estimate(size_t exsig, BasisData bd) {
    return nrow_estimate(decode_nfrm_cre(exsig), decode_nfrm_ann(exsig),
                         decode_nbos_cre(exsig), decode_nbos_ann(exsig), bd);
}

Rdm::Rdm(const fciqmc_config::Rdms &opts, size_t ranksig, BasisData bd, size_t nelec, size_t nvalue) :
        Communicator<MaeRow, MaeRow, true>(
                "rdm_" + to_string(ranksig), nrow_estimate(ranksig, bd),
                nrow_estimate(ranksig, bd), opts.m_buffers, opts.m_load_balancing,
                {{ranksig, nvalue}}, {{ranksig, nvalue}}
        ),
        m_ranksig(ranksig), m_rank(decode_nfrm_cre(ranksig)),
        m_nfrm_cre(decode_nfrm_cre(ranksig)), m_nfrm_ann(decode_nfrm_ann(ranksig)),
        m_nbos_cre(decode_nbos_cre(ranksig)), m_nbos_ann(decode_nbos_ann(ranksig)),
        m_lookup_inds(ranksig), m_bd(bd) {

    /*
     * if contributing exsig != ranksig, there is promotion to do
     * the promoter to use is given by the difference between either fermion element of the ranksig and that of the
     * contributing exsig
     */
    const auto rank = decode_nfrm_cre(m_ranksig);
    m_frm_promoters.reserve(rank + 1);
    for (size_t nins = 0ul; nins <= rank; ++nins)
        m_frm_promoters.emplace_back(nelec + nins - rank, nins);
}

void Rdm::make_contribs(const field::FrmOnv &src_onv, const conn::FrmOnv &conn,
                        const FrmOps &com, const wf_t &contrib) {
    const auto exlvl = conn.m_cre.size();
    DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                      "this method should not have been delegated given the exsig of the contribution");
    /*
     * number of "inserted" fermion creation/annihilation operator pairs
     */
    const auto nins = m_rank - exlvl;
    /*
     * this determines the precomputed promoter required
     */
    const auto &promoter = m_frm_promoters[nins];
    /*
     * apply each combination of the promoter deterministically
     */
    for (size_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
        auto phase = promoter.apply(icomb, conn, com, m_lookup_inds.m_frm);

        auto irank_send = m_ra.get_rank(m_lookup_inds);
        DEBUG_ASSERT_TRUE(m_lookup_inds.is_ordered(),
            "operators of each kind should be stored in ascending order of their orbital (or mode) index");
        auto &send_table = send(irank_send);
        size_t irow = *send_table[m_lookup_inds];
        if (irow == ~0ul) irow = send_table.insert(m_lookup_inds);
        send_table.m_row.jump(irow);
        /*
         * include the Fermi phase of the excitation
         */
        phase = phase ^ conn.phase(src_onv);
        send_table.m_row.m_values[0] += phase ? -contrib : contrib;
    }
}

void Rdm::make_contribs(const field::FrmBosOnv &src_onv, const conn::FrmBosOnv &conn,
                        const com_ops::FrmBos &com, const wf_t &contrib) {
    auto exsig = conn.exsig();
    m_lookup_inds.zero();
    if (is_pure_frm(exsig) && is_pure_frm(m_ranksig))
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib);
    /*
     * fermion promotion (if any) is handled in the delegated method, but if this is a hopping-coupled or density-coupled
     * boson contribution, then the boson occupation factor must also be included (like we have to consider in the
     * matrix elements in FrmBosHamiltonian)
     */
    if (decode_nbos(exsig)==1) {
        // "ladder" contribution
        m_lookup_inds = conn.m_bos;
        auto occ_fac = conn.m_bos.occ_fac(src_onv.m_bos);
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib * occ_fac);
    }

    m_lookup_inds.m_frm.zero();
    if (m_nbos_cre == 1ul && m_nbos_ann == 1ul && is_pure_frm(exsig)) {
        /*
         * this is the only currently supported situation in which boson promotion is required: a purely fermionic
         * (nbos_cre = 0, nbos_ann = 0) excitation or a diagonal (!exsig) contribution
         * loop through all modes in src_onv to extract all "common" modes
         */
        for (size_t imode = 0ul; imode < src_onv.m_bos.nelement(); ++imode) {
            const auto ncom = src_onv.m_bos[imode];
            if (!ncom) continue;
            m_lookup_inds.m_bos.m_cre[0] = imode;
            m_lookup_inds.m_bos.m_ann[0] = imode;
            make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, double(ncom) * contrib);
        }
    }
}

void Rdm::end_cycle() {
    if (!send().buffer_size()) return;
    communicate();
    auto &row = m_comm.recv().m_row;
    if (!m_comm.recv().m_hwm) return;
    for (row.restart(); row.in_range(); row.step()) {
        auto irow_store = *m_store[row.m_inds];
        if (irow_store == ~0ul) irow_store = m_store.insert(row.m_inds);
        m_store.m_row.jump(irow_store);
        m_store.m_row.m_values += row.m_values;
    }
    m_comm.recv().clear();
}

void Rdm::save(hdf5::GroupWriter &gw) const {
    m_store.save(gw, to_string(m_ranksig));
}

std::array<defs::inds, defs::nexsig> Rdms::make_exsig_ranks() const {
    std::array<defs::inds, defs::nexsig> exsig_ranks;
    for (const auto &ranksig: m_active_ranksigs) {
        auto nfrm_cre = decode_nfrm_cre(ranksig);
        auto nfrm_ann = decode_nfrm_cre(ranksig);
        while (nfrm_cre != ~0ul && nfrm_ann != ~0ul) {
            auto nbos_cre = decode_nbos_cre(ranksig);
            auto nbos_ann = decode_nbos_ann(ranksig);
            while (nbos_cre != ~0ul && nbos_ann != ~0ul) {
                auto exsig = encode(nfrm_cre, nfrm_ann, nbos_cre, nbos_ann);
                DEBUG_ASSERT_LT(exsig, defs::nexsig, "exsig OOB");
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

Rdms::Rdms(const fciqmc_config::Rdms &opts, defs::inds ranksigs, BasisData bd, size_t nelec, const Epoch &accum_epoch) :
        Archivable("rdms", opts.m_archivable),
        m_active_ranksigs(std::move(ranksigs)), m_exsig_ranks(make_exsig_ranks()),
        m_work_conns(bd), m_work_com_ops(bd), m_explicit_ref_conns(opts.m_explicit_ref_conns),
        m_accum_epoch(accum_epoch) {
    for (const auto &ranksig: m_active_ranksigs) {
        REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
        REQUIRE_TRUE(conserves_nfrm(ranksig), "fermion non-conserving RDMs are not yet supported");
        REQUIRE_LE(decode_nbos_cre(ranksig), 1ul,
                   "RDMs with more than one boson creation operator are not yet supported");
        REQUIRE_LE(decode_nbos_ann(ranksig), 1ul,
                   "RDMs with more than one boson annihilation operator are not yet supported");
        REQUIRE_TRUE(m_rdms[ranksig] == nullptr, "No RDM rank should appear more than once in the specification");
        m_rdms[ranksig] = mem_utils::make_unique<Rdm>(opts, ranksig, bd, nelec, 1ul);
    }
    m_total_norm.m_local = 0.0;
}

Rdms::operator bool() const {
    return !m_active_ranksigs.empty();
}

bool Rdms::takes_contribs_from(const size_t &exsig) const {
    if (exsig > defs::nexsig) return false;
    return !m_exsig_ranks[exsig].empty();
}

void Rdms::make_contribs(const Mbf &src_onv, const conn::Mbf &conn, const com_ops::Mbf &com, const wf_t &contrib) {
    auto exsig = conn.exsig();
    if (!exsig) m_total_norm.m_local+=contrib;
    for (auto ranksig: m_exsig_ranks[exsig]) m_rdms[ranksig]->make_contribs(src_onv, conn, com, contrib);
}

void Rdms::make_contribs(const Mbf &src_onv, const Mbf &dst_onv, const wf_t &contrib) {
    m_work_conns[src_onv].connect(src_onv, dst_onv, m_work_com_ops[src_onv]);
    make_contribs(src_onv, m_work_conns[src_onv], m_work_com_ops[src_onv], contrib);
}

void Rdms::make_contribs(const SpawnTableRow &recv_row, const WalkerTableRow &dst_row, const Propagator &prop) {
    DEBUG_ASSERT_EQ(recv_row.m_dst_mbf, dst_row.m_mbf, "found row doesn't correspond to spawned dst");
    auto ipart_replica = dst_row.ipart_replica(recv_row.m_ipart_dst);
    defs::wf_t contrib = dst_row.m_weight[ipart_replica];
    // recover pre-death value of replica population (on average)
    contrib /= 1.0 - prop.tau() * (dst_row.m_hdiag - prop.m_shift.m_values[ipart_replica]);
    contrib = consts::conj(contrib);
    contrib *= recv_row.m_src_weight;
    make_contribs(recv_row.m_src_mbf, dst_row.m_mbf, contrib);
}

bool Rdms::all_stores_empty() const {
    for (auto &ranksig: m_active_ranksigs)
        if (!m_rdms[ranksig]->m_store.is_cleared())
            return false;
    return true;
}

void Rdms::end_cycle() {
    for (auto &ranksig: m_active_ranksigs) m_rdms[ranksig]->end_cycle();
    m_total_norm.all_sum();
}

bool Rdms::is_energy_sufficient(const Hamiltonian &ham) const {
    if (ham.m_bos.m_contribs_0011.is_nonzero(0ul)){
        if (!m_rdms[exsig_utils::ex_0011]) return false;
    }
    if (ham.m_frmbos.m_contribs_0010.is_nonzero(exsig_utils::ex_0010)){
        if (!m_rdms[exsig_utils::ex_0010]) return false;
    }
    if (ham.m_frmbos.m_contribs_0001.is_nonzero(exsig_utils::ex_0001)){
        if (!m_rdms[exsig_utils::ex_0001]) return false;
    }
    if (!m_rdms[exsig_utils::ex_double]) return false;
    return true;
}


defs::ham_comp_t Rdms::get_energy(const FrmHam& ham) const {
    if (ham.disabled()) return 0.0;
    auto& rdm2 = m_rdms[ex_2200];
    REQUIRE_TRUE_ALL(rdm2!=nullptr, "cannot compute energy without the 2RDM");
    defs::ham_t e1 = 0.0;
    defs::ham_t e2 = 0.0;
    defs::wf_t trace = 0.0;
    auto& row = rdm2->m_store.m_row;

    for (row.restart(); row.in_range(); row.step()){
        const size_t i=row.m_inds.m_frm.m_cre[0];
        const size_t j=row.m_inds.m_frm.m_cre[1];
        DEBUG_ASSERT_LT(i, j, "spin orbital creation indices should be ordered");
        const size_t k=row.m_inds.m_frm.m_ann[0];
        const size_t l=row.m_inds.m_frm.m_ann[1];
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
    // scale the one-body contribution by the number of two-body contributions
    e1 /= ham.m_nelec-1;
    e1 = mpi::all_sum(e1);
    e2 = mpi::all_sum(e2);
    trace = mpi::all_sum(trace);
    DEBUG_ASSERT_GT(std::abs(trace), 1e-14, "RDM trace should be non-zero");
    const auto norm = consts::real(trace) / integer_utils::nspair(ham.m_nelec);
    REQUIRE_NEARLY_EQ(norm / m_total_norm.m_reduced, 1.0, 1e-8,
                 "2RDM norm should match total of sampled diagonal contributions");
    return consts::real(ham.m_e_core) + (consts::real(e1) + consts::real(e2)) / norm;
}

defs::ham_comp_t Rdms::get_energy(const FrmBosHam &ham, size_t nelec, size_t exsig) const {
    if (ham.disabled()) return 0.0;
    auto& rdm = m_rdms[exsig];
    REQUIRE_TRUE_ALL(exsig_utils::decode_nbos(exsig)==1,
                     "currently only supporting linear boson operators in the ladder term");
    REQUIRE_TRUE_ALL(rdm!=nullptr, "cannot compute energy without the "+exsig_utils::to_string(exsig)+"-RDM");
    defs::ham_t e_uncoupled = 0.0; // 0001 and 0010
    defs::ham_t e_coupled = 0.0; // 1101 and 1110
    auto& row = rdm->m_store.m_row;
    bool cre = exsig_utils::decode_nbos_cre(exsig);

    for (row.restart(); row.in_range(); row.step()){
        const size_t p=row.m_inds.m_frm.m_cre[0];
        const size_t q=row.m_inds.m_frm.m_ann[0];
        const size_t n=cre ? row.m_inds.m_bos.m_cre[0] : row.m_inds.m_bos.m_ann[0];

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
    REQUIRE_NEARLY_EQ(consts::imag(e), 0.0, 1e-12, "energy should be purely real");
    return consts::real(e);
}

defs::ham_comp_t Rdms::get_energy(const BosHam &ham) const {
    if (ham.disabled()) return 0.0;
    auto& rdm = m_rdms[exsig_utils::ex_0011];
    REQUIRE_TRUE_ALL(rdm!=nullptr, "cannot compute energy without the 0011-RDM");
    defs::ham_t e = 0.0;
    auto& row = rdm->m_store.m_row;

    for (row.restart(); row.in_range(); row.step()){
        const size_t n=row.m_inds.m_bos.m_cre[0];
        const size_t m=row.m_inds.m_bos.m_ann[0];
        REQUIRE_EQ(n, m, "0011-RDM should currently only take 0000-exsig contributions");
        const auto rdm_element = row.m_values[0];
        e += rdm_element*ham.get_coeff_0011(n, m);
    }
    e = mpi::all_sum(e) / m_total_norm.m_reduced;
    REQUIRE_TRUE(consts::nearly_zero(consts::imag(e), 1e-12), "energy should be purely real")
    return consts::real(e);
}
