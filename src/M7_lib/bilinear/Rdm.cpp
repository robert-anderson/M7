//
// Created by Robert J. Anderson on 11/08/2021.
//

#include "Rdm.h"
#include "M7_lib/util/SmartPtr.h"
#include "SpinfreeRdm.h"

uint_t Rdm::nrow_estimate(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann, sys::Size basis_size) {
    double nrow = 1.0;
    nrow *= integer::combinatorial(basis_size.m_frm.m_nspinorb, nfrm_cre);
    nrow *= integer::combinatorial(basis_size.m_frm.m_nspinorb, nfrm_ann);
    nrow *= integer::combinatorial(basis_size.m_bos, nbos_cre);
    nrow *= integer::combinatorial(basis_size.m_bos, nbos_ann);
    nrow /= integer::factorial(nfrm_cre + nfrm_ann);
    nrow /= integer::factorial(nbos_cre + nbos_ann);
    return nrow;
}

uint_t Rdm::nrow_estimate(uint_t exsig, sys::Size basis_size) {
    return nrow_estimate(decode_nfrm_cre(exsig), decode_nfrm_ann(exsig),
                         decode_nbos_cre(exsig), decode_nbos_ann(exsig), basis_size);
}

str_t Rdm::name(str_t str, uint_t ranksig) const {
    return str.empty() ? exsig::to_string(ranksig) : str;
}

void Rdm::add_to_send_table(const MaeInds &inds, wf_t contrib) {
    DEBUG_ASSERT_EQ(inds.m_exsig, m_indsig, "incorrect number of operators in MAE indices");
    const auto irank_send = m_dist.irank(inds);
    DEBUG_ASSERT_TRUE(inds.is_ordered(),
                      "operators of each kind should be stored in ascending order of their orbital (or mode) index");
    auto &send_table = send(irank_send);
    send_table.associate(m_send_row);
    auto lookup = send_table.lookup(inds, m_send_row);
    if (!lookup) send_table.insert(inds, m_send_row);
    m_send_row.m_values[0] += contrib;
}

void Rdm::frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                            const FrmOps& com, const wf_t& contrib) {
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
    const auto& promoter = m_frm_promoters[nins];
    /*
     * apply each combination of the promoter deterministically
     */
    for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
        auto phase = promoter.apply(icomb, conn, com, m_full_inds.m_frm);
        /*
         * include the Fermi phase of the excitation
         */
        phase ^= conn.phase(src_onv);
        add_to_send_table(m_full_inds, phase ? -contrib : contrib);
    }
}

void Rdm::frmbos_make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                               const com_ops::FrmBos& com, const wf_t& contrib) {
    auto exsig = conn.exsig();
    m_full_inds.zero();
    if (is_pure_frm(exsig) && is_pure_frm(m_ranksig))
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib);
    /*
     * fermion promotion (if any) is handled in the delegated method, but if this is a hopping-coupled or density-coupled
     * boson contribution, then the boson occupation factor must also be included (like we have to consider in the
     * matrix elements in FrmBosHamiltonian)
     */
    if (decode_nbos(exsig)==1) {
        // "ladder" contribution
        m_full_inds = conn.m_bos;
        auto occ_fac = src_onv.m_bos.occ_fac(conn.m_bos);
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib * occ_fac);
    }

    m_full_inds.m_frm.zero();
    if (m_nbos_cre == 1ul && m_nbos_ann == 1ul && is_pure_frm(exsig)) {
        /*
         * this is the only currently supported situation in which boson promotion is required: a purely fermionic
         * (nbos_cre = 0, nbos_ann = 0) excitation or a diagonal (!exsig) contribution
         * loop through all modes in src_onv to extract all "common" modes
         */
        for (uint_t imode = 0ul; imode < src_onv.m_bos.nelement(); ++imode) {
            const auto ncom = src_onv.m_bos[imode];
            if (!ncom) continue;
            m_full_inds.m_bos.m_cre[0] = imode;
            m_full_inds.m_bos.m_ann[0] = imode;
            make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, double(ncom) * contrib);
        }
    }
}

Rdm::Rdm(uint_t ranksig, uint_t indsig, sys::Size basis_size, uint_t nelec, uint_t nvalue,
         DistribOptions dist_opts, Sizing store_sizing, Sizing comm_sizing, str_t name) :
        communicator::MappedSend<MaeRow, MaeRow>(
                "rdm_" + (name.empty() ? this->name(name, ranksig) : name),
                // distributed "store" table row
                MaeRow(indsig, nvalue),
                dist_opts, store_sizing,
                // send/recv table row
                MaeRow(indsig, nvalue),
                comm_sizing),
        m_basis_size(basis_size), m_ranksig(ranksig), m_rank(decode_nfrm_cre(ranksig)),
        m_nfrm_cre(decode_nfrm_cre(ranksig)), m_nfrm_ann(decode_nfrm_ann(ranksig)),
        m_nbos_cre(decode_nbos_cre(ranksig)), m_nbos_ann(decode_nbos_ann(ranksig)),
        m_indsig(indsig), m_rank_ind(decode_nfrm_cre(m_indsig)),
        m_nfrm_cre_ind(decode_nfrm_cre(m_indsig)), m_nfrm_ann_ind(decode_nfrm_ann(m_indsig)),
        m_nbos_cre_ind(decode_nbos_cre(m_indsig)), m_nbos_ann_ind(decode_nbos_ann(m_indsig)),
        m_full_inds(ranksig), m_uncontracted_inds(m_indsig), m_name(name),
        m_send_row(m_send_recv.m_row), m_recv_row(m_send_recv.m_row), m_store_row(m_store.m_row), m_nelec(nelec) {
    /*
     * if contributing exsig != ranksig, there is promotion to do
     * the promoter to use is given by the difference between either fermion element of the ranksig and that of the
     * contributing exsig
     */
    m_frm_promoters.reserve(m_rank + 1);
    for (uint_t nins = 0ul; nins <= m_rank; ++nins)
        m_frm_promoters.emplace_back(nelec + nins - m_rank, nins);
}


void Rdm::end_cycle() {
    if (!send().buffer_size()) return;
    communicate();
    if (!m_send_recv.recv().m_hwm) return;
    for (m_recv_row.restart(); m_recv_row.in_range(); m_recv_row.step()) {
        auto lookup = m_store.lookup(m_recv_row.m_inds, m_store_row);
        if (!lookup) m_store.insert(m_recv_row.m_inds, m_store_row);
        m_store_row.m_values += m_recv_row.m_values;
    }
    m_send_recv.recv().clear();
}

void Rdm::save(hdf5::NodeWriter& gw) const {
    m_store.save(gw, this->name());
}

FockRdm4::FockRdm4(const conf::Rdms &opts, sys::Size basis_size, uint_t nelec, uint_t nvalue) :
        Rdm(opts, exsig::ex_4400, exsig::ex_3300, basis_size, nelec, nvalue, "4400f"),
        m_fock(basis_size.m_frm.m_nsite){
    logging::info("loading generalized Fock matrix for CASPT2 contracted 4RDM accumulation");
    hdf5::FileReader reader(opts.m_fock_4rdm.m_fock_path);
    REQUIRE_TRUE(reader.child_exists("ACT_FOCK_INDEX"), "invalid fock matrix file contents");
    REQUIRE_TRUE(reader.child_exists("ACT_FOCK_VALUES"), "invalid fock matrix file contents");
    v_t<int64_t> inds;
    reader.read_data("ACT_FOCK_INDEX", inds);
    v_t<double> values;
    reader.read_data("ACT_FOCK_VALUES", values);
    REQUIRE_EQ(inds.size(), values.size()*2, "incorrect number of indices");
    for (uint_t i = 0ul; i<values.size(); ++i) {
        // offset by one to account for the Fortran indices in the HDF5 format
        m_fock(inds[2*i]-1, inds[2*i+1]-1) = values[i];
    }
}

void FockRdm4::frm_make_contribs(const FrmOnv &src_onv, const conn::FrmOnv &conn,
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
    const auto& promoter = m_frm_promoters[nins];
    /*
     * apply each combination of the promoter deterministically
     */
    for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
        auto phase = promoter.apply(icomb, conn, com, m_full_inds.m_frm);
        /*
         * include the Fermi phase of the excitation
         */
        phase = phase ^ conn.phase(src_onv);
        bool contract_phase = true;
        // TODO: diagonal Fock optimisation
        for (uint_t icre_contract = 0ul; icre_contract < 4ul; ++icre_contract){
            const uint_t icre = m_full_inds.m_frm.m_cre[icre_contract];
            const auto isite_cre = src_onv.m_basis.isite(icre);
            const auto ispin_cre = src_onv.m_basis.ispin(icre);
            for (uint_t iann_contract = 0ul; iann_contract < 4ul; ++iann_contract) {
                const uint_t iann = m_full_inds.m_frm.m_ann[iann_contract];
                const auto isite_ann = src_onv.m_basis.isite(iann);
                const auto ispin_ann = src_onv.m_basis.ispin(iann);

                if (ispin_ann != ispin_cre) continue;

                const auto fock_element = m_fock(isite_cre, isite_ann);

                contract_phase = !contract_phase;
                /*
                 * fill the uncontracted indices (those identifying the elements of the intermediate)
                 */
                for (uint_t iuncontract = 0ul; iuncontract < 4ul; ++iuncontract) {
                    // the next position in the uncontracted indices to fill
                    uint_t i;
                    if (iuncontract!=icre_contract) {
                        i = iuncontract - (iuncontract > icre_contract ? 1 : 0);
                        m_uncontracted_inds.m_frm.m_cre[i] = m_full_inds.m_frm.m_cre[iuncontract];
                    }
                    if (iuncontract!=iann_contract) {
                        i = iuncontract - (iuncontract > iann_contract ? 1 : 0);
                        m_uncontracted_inds.m_frm.m_ann[i] = m_full_inds.m_frm.m_ann[iuncontract];
                    }
                }

                /*
                 * include the Fermi phase of the contraction and the promotion when computing the overall sign
                 */
                add_to_send_table(m_uncontracted_inds,
                                  ((phase ^ contract_phase) ? -contrib : contrib) * fock_element);
            }
        }
    }
}


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

Rdms::Rdms(const conf::Rdms& opts, uintv_t ranksigs, sys::Size basis_size, uint_t nelec, const Epoch& accum_epoch) :
        Archivable("rdms", opts.m_archivable),
        m_spinfree(opts.m_spinfree), m_rdm_ranksigs(ranksigs), m_exsig_ranks(make_exsig_ranks()),
        m_work_conns(basis_size), m_work_com_ops(basis_size), m_explicit_ref_conns(opts.m_explicit_ref_conns),
        m_accum_epoch(accum_epoch), m_nelec(nelec) {
    for (const auto& ranksig: ranksigs) {
        REQUIRE_TRUE(ranksig, "multidimensional estimators require a nonzero number of SQ operator indices");
        REQUIRE_TRUE(conserves_nfrm(ranksig), "fermion non-conserving RDMs are not yet supported");
        REQUIRE_LE(decode_nbos_cre(ranksig), 1ul,
                   "RDMs with more than one boson creation operator are not yet supported");
        REQUIRE_LE(decode_nbos_ann(ranksig), 1ul,
                   "RDMs with more than one boson annihilation operator are not yet supported");
        REQUIRE_TRUE(m_rdms[ranksig] == nullptr, "No RDM rank should appear more than once in the specification");
        m_rdms[ranksig] = smart_ptr::make_unique<Rdm>(opts, ranksig, ranksig, basis_size, nelec, 1ul);
    }
    if (opts.m_fock_4rdm.m_enabled) m_fock_rdm4 = smart_ptr::make_unique<FockRdm4>(opts, basis_size, nelec, 1ul);
    m_total_norm.m_local = 0.0;
}

Rdms::operator bool() const {
    return !m_rdm_ranksigs.empty() || m_fock_rdm4;
}

bool Rdms::takes_contribs_from(uint_t exsig) const {
    if (exsig > exsig::c_ndistinct) return false;
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

void Rdms::make_contribs(const Spawn& recv_row, const Walker& dst_row, const Propagator& prop) {
    DEBUG_ASSERT_EQ(recv_row.m_dst_mbf, dst_row.m_mbf, "found row doesn't correspond to spawned dst");
    auto ipart_replica = dst_row.ipart_replica(recv_row.m_ipart_dst);
    wf_t contrib = dst_row.m_weight[ipart_replica];
    // recover pre-death value of replica population (on average)
    contrib /= 1.0 - prop.tau() * (dst_row.m_hdiag - prop.m_shift.m_values[ipart_replica]);
    contrib = arith::conj(contrib);
    contrib *= recv_row.m_src_weight;
    make_contribs(recv_row.m_src_mbf, dst_row.m_mbf, contrib);
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

    for (row.restart(); row.in_range(); row.step()){
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
    // scale the one-body contribution by the number of two-body contributions
    e1 /= m_nelec-1;
    e1 = mpi::all_sum(e1);
    e2 = mpi::all_sum(e2);
    trace = mpi::all_sum(trace);
    DEBUG_ASSERT_GT(std::abs(trace), 1e-14, "RDM trace should be non-zero");
    const auto norm = arith::real(trace) / integer::nspair(m_nelec);
    REQUIRE_NEARLY_EQ(norm, m_total_norm.m_reduced, "2RDM norm should match total of sampled diagonal contributions");
    return arith::real(ham.m_e_core) + (arith::real(e1) + arith::real(e2)) / norm;
}

ham_comp_t Rdms::get_energy(const FrmBosHam& /*ham*/, uint_t /*nelec*/, uint_t /*exsig*/) const {
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

    for (row.restart(); row.in_range(); row.step()){
        const uint_t n=row.m_inds.m_bos.m_cre[0];
        const uint_t m=row.m_inds.m_bos.m_ann[0];
        REQUIRE_EQ(n, m, "0011-RDM should currently only take 0000-exsig contributions");
        const auto rdm_element = row.m_values[0];
        e += rdm_element*ham.get_coeff_0011(n, m);
    }
    e = mpi::all_sum(e) / m_total_norm.m_reduced;
    REQUIRE_TRUE(fptol::numeric_real(e), "energy should be purely real")
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
