//
// Created by Robert J. Anderson on 11/08/2021.
//

#include "Rdm.h"
#include "M7_lib/util/Pointer.h"
#include "SpinfreeRdm.h"

uint_t Rdm::nrow_estimate(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann, sys::Size basis_size) {
    double nrow = 1.0;
    nrow *= integer::combinatorial(basis_size.m_frm.m_nspinorb, nfrm_cre);
    nrow *= integer::combinatorial(basis_size.m_frm.m_nspinorb, nfrm_ann);
    nrow *= integer::combinatorial_with_repetition(basis_size.m_bos, nbos_cre);
    nrow *= integer::combinatorial_with_repetition(basis_size.m_bos, nbos_ann);
    nrow /= integer::factorial(nfrm_cre + nfrm_ann);
    nrow /= integer::factorial(nbos_cre + nbos_ann);
    return nrow;
}

uint_t Rdm::nrow_estimate(OpSig exsig, sys::Size basis_size) {
    return nrow_estimate(exsig.nfrm_cre(), exsig.nfrm_ann(), exsig.nbos_cre(), exsig.nbos_ann(), basis_size);
}

str_t Rdm::name(str_t str, OpSig ranksig) const {
    return str.empty() ? ranksig.to_string() : str;
}

void Rdm::add_to_send_table(const MaeInds &inds, wf_t contrib) {
    DEBUG_ASSERT_EQ(inds.m_exsig, m_indsig, "incorrect number of operators in MAE indices");
    const auto irank_send = m_dist.irank(inds);
    DEBUG_ASSERT_TRUE(!m_ordered_inds || inds.is_ordered(),
                      "operators of each kind should be stored in ascending order of their orbital (or mode) index");
    auto &send_table = send(irank_send);
    send_table.associate(m_send_row);
    auto lookup = send_table.lookup(inds, m_send_row);
    if (!lookup) send_table.insert(inds, m_send_row);
    m_send_row.m_values[0] += contrib;
}

Rdm::Rdm(OpSig ranksig, OpSig indsig, sys::Sector sector, uint_t nvalue,
         DistribOptions dist_opts, Sizing store_sizing, Sizing comm_sizing, str_t name) :
        communicator::MappedSend<MaeRow, MaeRow>(
                "rdm_" + (name.empty() ? this->name(name, ranksig) : name),
                // distributed "store" table row
                MaeRow(indsig, nvalue),
                dist_opts, store_sizing,
                // send/recv table row
                MaeRow(indsig, nvalue),
                comm_sizing),
        m_sector(sector), m_ranksig(ranksig),
        m_nfrm_cre(ranksig.nfrm_cre()), m_nfrm_ann(ranksig.nfrm_ann()),
        m_nbos_cre(ranksig.nbos_cre()), m_nbos_ann(ranksig.nbos_ann()),
        m_indsig(indsig), //m_rank_ind(m_indsig.(m_indsig)),
        m_nfrm_cre_ind(m_indsig.nfrm_cre()), m_nfrm_ann_ind(m_indsig.nfrm_ann()),
        m_nbos_cre_ind(m_indsig.nbos_cre()), m_nbos_ann_ind(m_indsig.nbos_ann()),
        m_full_inds(ranksig), m_uncontracted_inds(m_indsig), m_name(name),
        m_send_row(m_send_recv.m_row), m_recv_row(m_send_recv.m_row), m_store_row(m_store.m_row) {
    /*
     * if contributing exsig != ranksig, there is promotion to do
     * the promoter to use is given by the difference between either fermion element of the ranksig and that of the
     * contributing exsig
     */
    const auto frm_rank = m_ranksig.nfrm_cre();
    REQUIRE_TRUE(m_ranksig.conserves_nfrm(), "RDMs not implemented for non-conserved fermion number");
    m_frm_promoters.reserve(frm_rank + 1);
    for (uint_t nins = 0ul; nins <= frm_rank; ++nins) {
        const auto nexcit = frm_rank - nins;
        m_frm_promoters.emplace_back(sector.m_frm.m_elecs - nexcit, opsig::frm(nexcit), nins);
    }
}

Rdm::Rdm(const conf::Rdms& opts, OpSig ranksig, OpSig indsig, sys::Sector sector, uint_t nvalue, str_t name) :
        Rdm(ranksig, indsig, sector, nvalue, opts.m_distribution,
            // store sizing
            Sizing{nrow_estimate(indsig, sector.size()), opts.m_buffers.m_store_exp_fac},
            // send/recv sizing
            Sizing{nrow_estimate(indsig, sector.size()), opts.m_buffers.m_comm_exp_fac}, name){}

void Rdm::end_cycle() {
    if (!send().buffer_size()) return;
    communicate();
    if (m_send_recv.recv().empty()) return;
    for (m_recv_row.restart(); m_recv_row; ++m_recv_row) {
        auto lookup = m_store.lookup(m_recv_row.m_inds, m_store_row);
        if (!lookup) m_store.insert(m_recv_row.m_inds, m_store_row);
        m_store_row.m_values += m_recv_row.m_values;
    }
    m_send_recv.recv().clear();
}

void Rdm::save(hdf5::NodeWriter& gw) const {
    m_store.save(gw, this->name(), true);
}

void PureRdm::frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn, const FrmOps& com, wf_t contrib) {
    const auto exlvl = conn.m_cre.size();
    DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                      "this method should not have been delegated given the exsig of the contribution");
    const auto rank = m_ranksig.nfrm_cre();
    /*
     * number of "inserted" fermion creation/annihilation operator pairs
     */
    const auto nins = rank - exlvl;
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

void PureRdm::frmbos_make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                               const com_ops::FrmBos& com, wf_t contrib) {
    auto exsig = conn.exsig();
    m_full_inds.zero();
    if (exsig.is_pure_frm() && m_ranksig.is_pure_frm()) make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib);
    /*
     * fermion promotion (if any) is handled in the delegated method, but if this is a hopping-coupled or density-coupled
     * boson contribution, then the boson occupation factor must also be included (like we have to consider in the
     * matrix elements in FrmBosHamiltonian)
     */
    if (exsig.nbos()==1) {
        // "ladder" contribution
        m_full_inds = conn.m_bos;
        auto occ_fac = src_onv.m_bos.occ_fac(conn.m_bos);
        make_contribs(src_onv.m_frm, conn.m_frm, com.m_frm, contrib * occ_fac);
    }

    m_full_inds.m_frm.zero();
    if (m_nbos_cre == 1ul && m_nbos_ann == 1ul && exsig.is_pure_frm()) {
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
