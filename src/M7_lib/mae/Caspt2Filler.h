//
// Created by Robert John Anderson on 22/08/2023.
//

#ifndef M7_CASPT2FILLER_H
#define M7_CASPT2FILLER_H

#include "M7_lib/wavefunction/WalkerTable.h"
#include "M7_lib/communication/Communicator.h"
#include "M7_lib/field/Mbf.h"
#include "M7_lib/bilinear/Rdm.h"
#include "M7_lib/bilinear/SpinfreeRdm.h"
#include "M7_lib/bilinear/FockRdm4.h"
#include "M7_lib/bilinear/Rdms.h"
#include "M7_lib/table/Smuvi.h"

#if 0
namespace pose {
    struct Row : public ::Row {
        field::Numbers<uint8_t, 1> m_inds;
        field::Numbers<wf_t, 1> m_values;
        field::Numbers<uint8_t, 1> &key_field() {
            return m_inds;
        }

        Row() :
            m_inds(this, {4ul}, "indices"),
            m_values(this, {1ul}, "values"){}
    };

    class Rdm : public communicator::MappedSend<Row, Row> {

        conn::Mbf m_work_conn;

        buffered::Numbers<uint8_t, 1> m_work_inds;

        void add_to_send_table(const field::Numbers<uint8_t, 1> &inds, wf_t contrib) {
            const auto irank = m_dist.irank(inds);
            auto& send_table = send(irank);
            auto& row = send_table.lookup_or_insert(inds);
            row.m_values += contrib;
        }

        void add_to_send_table(uintp_t p1, uintp_t p2, wf_t contrib) {
            m_work_inds[0] = p1.first;
            m_work_inds[1] = p1.second;
            m_work_inds[2] = p2.first;
            m_work_inds[3] = p2.second;
            add_to_send_table(m_work_inds, contrib);
        }

        void make_contribs(const field::Mbf& mbf, wf_t contrib) {
            /*
             * I is identified with the mbf argument.
             * let ~I be the set of spinorbs not occupied I.
             * make contribs of the form:
             * 1. <I| e_ii e_jj |I> for all i and j in I
             * 2. <I| e_ij e_ji |I> for all i in I, and all j in ~I where s(i) == s(j)
             */
            const auto& basis = mbf.m_basis;

            auto fill_1 = [&](uint_t i) {
                const auto oi = basis.isite(i);
                auto inner_fn = [&](uint_t j) {
                    const auto oj = basis.isite(j);
                    add_to_send_table({oi, oi}, {oj, oj}, contrib);
                };
                mbf.foreach_setbit(inner_fn);
            };
            mbf.foreach_setbit(fill_1);

            auto fill_2 = [&](uint_t i) {
                const auto oi = basis.isite(i);
                const auto si = basis.ispin(i);
                auto inner_fn = [&](uint_t j) {
                    const auto oj = basis.isite(j);
                    add_to_send_table({oj, oi}, {oi, oj}, contrib);
                };
                si ? mbf.foreach_beta(inner_fn) : mbf.foreach_alpha(inner_fn);
            };
            mbf.foreach_clrbit(fill_2);
        }

        void make_contribs(const field::Mbf& src, const field::Mbf& dst, uintp_t p1, wf_t contrib) {
            /*
             * let I and J be determinants which are the same apart from two spinorbs with the same spin s
             * i and j are the first and second elements of p1
             * i is in I but not J, j is in J but not I.
             *
             * I is identified with the dst argument.
             * J is identified with the src argument.
             *
             * let K be the set of spinorbs occupied in both I and J
             * let ~K be the set of spinorbs occupied in neither I nor J
             *
             * we certainly have contribs to the spin-resolved POSE 2RDM of these forms:
             * 1. <I| e_kk e_ij |J> for all k in I
             * 2. <I| e_ij e_kk |J> for all k in J
             *
             * we also have contribs of the form:
             * 3. <I| e_kj e_ik |J> for all k of spin s in K
             * 4. <I| e_ik e_kj |J> for all k of spin s in ~K
             *
             * the Fermi phases of these contribs relative to the normal-ordered e_ij connection are all +1.
             */
            const auto& i = p1.first;
            const auto& j = p1.second;
            DEBUG_ASSERT_TRUE(dst.get(i), "i should be occupied in the dst");
            DEBUG_ASSERT_FALSE(dst.get(j), "j should not be occupied in the dst");
            DEBUG_ASSERT_FALSE(src.get(i), "i should not be occupied in the src");
            DEBUG_ASSERT_TRUE(src.get(j), "j should be occupied in the src");

            const auto& basis = src.m_basis;
            const auto io = basis.isite(i);
            const auto jo = basis.isite(j);
            const auto spin = basis.ispin(io);
            if (basis.ispin(jo) != spin) return;

            auto fill_1 = [&](uint_t k) {
                const auto ko = basis.isite(k);
                add_to_send_table({ko, ko}, {io, jo}, contrib);
            };
            spin ? dst.foreach_beta(fill_1) : dst.foreach_alpha(fill_1);

            auto fill_2 = [&](uint_t k) {
                const auto ko = basis.isite(k);
                add_to_send_table({io, jo}, {ko, ko}, contrib);
            };
            spin ? src.foreach_beta(fill_2) : src.foreach_alpha(fill_2);

            auto fill_3 = [&](uint_t k) {
                if (k == i) return;
                const auto ko = basis.isite(k);
                add_to_send_table({ko, jo}, {io, ko}, contrib);
            };
            dst.foreach_setbit(fill_3);

            auto fill_4 = [&](uint_t k) {
                if (k == j) return;
                const auto ko = basis.isite(k);
                add_to_send_table({io, ko}, {ko, jo}, contrib);
            };
            dst.foreach_clrbit(fill_4);
        }

        void make_contribs(const field::Mbf& src, uintp_t p1, uintp_t p2, wf_t contrib) {
            const auto& basis = src.m_basis;
            auto s1 = basis.ispin(p1.first);
            auto s2 = basis.ispin(p2.first);
            DEBUG_ASSERT_EQ(s1, basis.ispin(p1.second), "first pair is not spin conserving");
            DEBUG_ASSERT_EQ(s2, basis.ispin(p2.second), "second pair is not spin conserving");
            /*
             * i, j = p1
             * k, l = p2
             *
             * 1. <I| e_ij e_kl |J>
             * 2. <I| e_kl e_ij |J> (+)
             * 3. <I| e_il e_kj |J> (-)
             * 4. <I| e_kj e_il |J> (-)
             *
             * 3 and 4 only contribute if s1==s2
             */
            const auto io = basis.isite(p1.first);
            const auto jo = basis.isite(p1.second);
            const auto ko = basis.isite(p2.first);
            const auto lo = basis.isite(p2.second);
            add_to_send_table({io, jo}, {ko, lo}, contrib);
            add_to_send_table({ko, lo}, {io, jo}, contrib);
            if (s1!=s2) return;
            add_to_send_table({io, lo}, {ko, jo}, contrib);
            add_to_send_table({ko, jo}, {io, lo}, contrib);
        }

    public:

        Rdm(str_t name, sys::Basis basis):
            communicator::MappedSend<Row, Row>(name, {}, DistribOptions(), Sizing{1000, 1.0}, {}, Sizing{1000ul, 1.0}),
            m_work_conn(basis.size()), m_work_inds({4ul}){}

        void make_contribs(const field::Mbf& src, const field::Mbf& dst, wf_t contrib) {
            const auto exsig = mbf::exsig(src, dst);
            if (exsig.nfrm_cre() > 2) return;
            const auto& basis = src.m_basis;
            m_work_conn.connect(src, dst);
            switch (exsig.nfrm_cre()) {
                case 0ul:
                    make_contribs(src, contrib);
                    break;
                case 1ul: {
                    if (m_work_conn.phase(src)) contrib*=-1.0;
                    make_contribs(src, dst, {m_work_conn.m_cre[0], m_work_conn.m_ann[0]}, contrib);
                    break;
                }
                case 2ul: {
                    const auto &i = m_work_conn.m_cre[0];
                    const auto &j = m_work_conn.m_cre[1];
                    const auto &k = m_work_conn.m_ann[0];
                    const auto &l = m_work_conn.m_ann[1];
                    if (m_work_conn.phase(src)) contrib*=-1.0;
                    /*
                     * phase refers to a normal-ordered connection:
                     * i^+ j^+ l k
                     * i^+ k j^+ l (+)
                     * i^+ l j^+ k (-)
                     */
                    basis.ispin(i)==basis.ispin(k) ?
                        make_contribs(src, {i, k}, {j, l}, contrib):
                        make_contribs(src, {i, l}, {j, k}, -contrib);
                }
            }
        }
    };
}
#endif


/**
 * CASPT2 requires 3RDM and 4RDM*Fock which are both 6-indexed tensors.
 * the approach here is based on a histogrammed WF |0> and a two-fold RDM filling procedure:
 *  1. contract |0> with the Fock operator to form |F>
 *  2. loop over E_pq |psi_0>, forming |pq>
 *
 * The 3RDM is computable from the transition 2RDM of |pq> and |0>,
 * and the 4RDM*Fock is computable from the transition 2RDM of |pq> and |F>
 *
 */
class Caspt2Filler {
    /**
     * Histogrammable set of determinants of which to compute the outer product in filling the MAEs
     */
    const Table<MbfWeightRow>& m_hist;
    /**
     * the result of F * m_hist where F = sum_pq f_pq E_pq
     */
    communicator::BasicSend<MbfWeightRow, MbfWeightRow> m_fock_x_hist;
    /**
     * normal ordered, spin-resolved RDMs being filled
     */
    Rdms* m_rdms = nullptr;
    /**
     * working objects for connections and common indices
     */
    mutable suite::Conns m_work_conns;
    mutable suite::ComOps m_work_com_ops;


    struct SpinChannelToIndsSmuvi : Smuvi<field::FrmOnvSpinChannel, field::Number<uint_t>> {
    private:
        buffered::Number<uint_t> m_inserter;

    public:
        SpinChannelToIndsSmuvi(str_t name, size_t nsite):
                Smuvi<field::FrmOnvSpinChannel, field::Number<uint_t>>(
                        std::move(name),
                        field::FrmOnvSpinChannel(nullptr, nsite),
                        field::Number<uint_t>(nullptr)){}

        void insert(const field::FrmOnvSpinChannel& key, uint_t value) {
            m_inserter = value;
            Smuvi<field::FrmOnvSpinChannel, field::Number<uint_t>>::insert(key, m_inserter);
        }
    };

    struct SpinChannelToSpinChannelSmuvi : Smuvi<field::FrmOnvSpinChannel, field::FrmOnvSpinChannel> {
        SpinChannelToSpinChannelSmuvi(str_t name, size_t nsite):
                Smuvi<field::FrmOnvSpinChannel, field::FrmOnvSpinChannel>(
                        std::move(name),
                        field::FrmOnvSpinChannel(nullptr, nsite),
                        field::FrmOnvSpinChannel(nullptr, nsite)){}
    };

    struct SpinChannelToFrmOnvSmuvi : Smuvi<field::FrmOnvSpinChannel, field::FrmOnv> {
        SpinChannelToFrmOnvSmuvi(str_t name, sys::frm::Basis basis):
                Smuvi<field::FrmOnvSpinChannel, field::FrmOnv>(
                        std::move(name),
                        field::FrmOnvSpinChannel(nullptr, basis.m_nsite),
                        field::FrmOnv(nullptr, basis)){}
    };
    /**
     *
     */
    SpinChannelToIndsSmuvi m_dets_contain_alpha;
    SpinChannelToIndsSmuvi m_dets_contain_beta;
    SpinChannelToSpinChannelSmuvi m_beta_with_alpha;
    SpinChannelToSpinChannelSmuvi m_alpha_with_beta;

    // Smuvi<field::FrmOnvSpinChannel> m_alpha_single_dict;
    // Smuvi<field::FrmOnvSpinChannel> m_beta_single_dict;

    // Smuvi<field::FrmOnvSpinChannel> m_alpha_singles;
    // Smuvi<field::FrmOnvSpinChannel> m_beta_singles;

    typedef std::pair<std::pair<uint_t, uint_t>, ham_t> pq_val_t;
    /**
     * make a linear combination psi1 = G * psi where G = sum_pq g_pq E_pq and psi is m_hist
     * @param psi1
     *  the output of G * psi
     * @param pq_vals
     *  a sparse representation of the elements of the contraction ((p, q), g_pq)
     */
    void make_psi1(communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const v_t<pq_val_t>& pq_vals) {
        psi1.m_store.clear();
        const auto displ = mpi::evenly_shared_displ(m_hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(m_hist.nrow_in_use());

        auto hist_row = m_hist.m_row;
        buffered::Mbf work_mbf(hist_row.m_mbf.m_basis);
        const auto& basis = work_mbf.m_basis;
        v_t<std::pair<conn::Mbf, ham_t>> conn_vals;
        for (auto& pq_val: pq_vals) {
            const auto p = pq_val.first.first;
            const auto q = pq_val.first.second;
            const auto contract_val = pq_val.second;
            for (uint_t ispin = 0ul; ispin < 2ul; ++ispin) {
                conn_vals.emplace_back(work_mbf, contract_val);
                if (p!=q) {
                    // not a diagonal
                    conn_vals.back().first.m_cre.set(basis.ispinorb(ispin, p));
                    conn_vals.back().first.m_ann.set(basis.ispinorb(ispin, q));
                }
            }
        }

        for (hist_row.restart(displ); hist_row.in_range(displ + count); ++hist_row) {
            for (auto& conn_val: conn_vals) {
                const auto& conn = conn_val.first;
                bool phase = false;
                Mbf* dst = &hist_row.m_mbf;
                if (conn.size()) {
                    if (mbf::destroys(conn, hist_row.m_mbf)) continue;
                    conn.apply(hist_row.m_mbf, work_mbf);
                    phase = conn.phase(hist_row.m_mbf);
                    dst = &work_mbf;
                }
                auto irank_dst = psi1.m_dist.irank(*dst);
                auto& send_row = psi1.m_send_recv.send(irank_dst).m_row;
                send_row.push_back_jump();
                send_row.m_mbf = *dst;
                send_row.m_weight = hist_row.m_weight;
                send_row.m_weight *= conn_val.second;
                if (phase) send_row.m_weight *= -1.0;
            }
        }

        psi1.communicate();
        /*
         * do mini annihilation loop
         */
        auto& recv_row = psi1.m_send_recv.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            auto& dst = psi1.m_store.lookup_or_insert(recv_row.m_mbf);
            dst.m_weight += recv_row.m_weight;
        }
    }

    /**
     * overload to make psi1 = G * psi where g_ij are elements of a dense matrix
     */
    void make_psi1(communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const dense::SquareMatrix<ham_t>& coeffs, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        auto fn = [&pq_vals](uint_t irow, uint_t icol, ham_t elem) {
            pq_vals.emplace_back(uintp_t(irow, icol), elem);
        };
        coeffs.foreach(fn, tol);
        make_psi1(psi1, pq_vals);
    }

    /**
     * overload to make psi1 = G * psi where g_ii are elements of a dense vector
     */
    void make_psi1(communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const dense::Vector<ham_t>& coeffs, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        for (uint_t ielem = 0ul; ielem < coeffs.nelement(); ++ielem)
            pq_vals.emplace_back(uintp_t(ielem, ielem), coeffs[ielem]);
        make_psi1(psi1, pq_vals);
    }

    // todo: make private again
public:

    void make_contribs(PureRdm* rdm, const field::Mbf& src, const field::Mbf& dst, wf_t contrib) const {
        auto& conn = m_work_conns[src];
        auto& com_ops = m_work_com_ops[src];
        conn.connect(src, dst, com_ops);
        rdm->make_contribs(src, conn, com_ops, contrib);
    }

    void fill_rdm1(PureRdm* rdm) const {
        REQUIRE_TRUE(rdm, "RDM pointer should not be null");
        REQUIRE_TRUE(rdm->m_ranksig == opsig::c_sing, "RDM object should be one-body");

        const auto displ = mpi::evenly_shared_displ(m_hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(m_hist.nrow_in_use());
        auto bra_row = m_hist.m_row;
        auto ket_row = bra_row;

        auto make_contrib_fn = [&]() {
            const auto contrib = bra_row.m_weight[0] * ket_row.m_weight[0];
            make_contribs(rdm, bra_row.m_mbf, ket_row.m_mbf, contrib);
        };

        buffered::FrmOnvSpinChannel alpha_channel(bra_row.m_mbf.m_basis.m_nsite);
        buffered::FrmOnvSpinChannel beta_channel(bra_row.m_mbf.m_basis.m_nsite);
        for (bra_row.restart(displ); bra_row.in_range(displ + count); ++bra_row) {
            bra_row.m_mbf.copy_alpha_to(alpha_channel);
            bra_row.m_mbf.copy_beta_to(beta_channel);
            // alpha-alpha
            m_dets_contain_beta.foreach_value(beta_channel, [&](const field::Number<uint_t>& iket){
                ket_row.jump(iket);
                const auto hamming_dist = bra_row.m_mbf.nalpha_not_in(ket_row.m_mbf);
                if (hamming_dist <= rdm->m_ranksig.nfrm_cre()) make_contrib_fn();
            });
            // beta-beta
            m_dets_contain_alpha.foreach_value(alpha_channel, [&](const field::Number<uint_t>& iket){
                ket_row.jump(iket);
                const auto hamming_dist = bra_row.m_mbf.nbeta_not_in(ket_row.m_mbf);
                if (hamming_dist <= rdm->m_ranksig.nfrm_cre() && hamming_dist > 0) make_contrib_fn();
            });
        }
    }

    void fill_rdm2(PureRdm* rdm) const {
        REQUIRE_TRUE(rdm, "RDM pointer should not be null");
        REQUIRE_TRUE(rdm->m_ranksig == opsig::c_doub, "RDM object should be two-body");

        const auto displ = mpi::evenly_shared_displ(m_hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(m_hist.nrow_in_use());
        auto bra_row = m_hist.m_row;
        auto ket_row = bra_row;

        auto make_contrib_fn = [&]() {
            const auto contrib = bra_row.m_weight[0] * ket_row.m_weight[0];
            make_contribs(rdm, bra_row.m_mbf, ket_row.m_mbf, contrib);
        };

        buffered::FrmOnvSpinChannel alpha_channel(bra_row.m_mbf.m_basis.m_nsite);
        buffered::FrmOnvSpinChannel beta_channel(bra_row.m_mbf.m_basis.m_nsite);
        for (bra_row.restart(displ); bra_row.in_range(displ + count); ++bra_row) {
            bra_row.m_mbf.copy_alpha_to(alpha_channel);
            bra_row.m_mbf.copy_beta_to(beta_channel);
        }

            // m_dets_contain_beta.foreach_entry_by_key(beta_channel, [&](uint_t iket){
            //     ket_row.jump(iket);
            //     const auto hamming_dist = bra_row.m_mbf.nalpha_not_in(ket_row.m_mbf);
            //     if (hamming_dist <= rdm->m_ranksig.nfrm_cre()) make_contrib_fn();
            // });

            // m_dets_contain_alpha.foreach_entry_by_key(alpha_channel, [&](uint_t iket){
            //     ket_row.jump(iket);
            //     const auto hamming_dist = bra_row.m_mbf.nbeta_not_in(ket_row.m_mbf);
            //     if (hamming_dist <= rdm->m_ranksig.nfrm_cre() && hamming_dist > 0) make_contrib_fn();
            // });

            /** Python
             *  for alpha_string in AlphaSingles[abra]:
             *     for beta_string in set(BetaWithAlpha[alpha_string]) & set(BetaSingles[bbra]):
             *         for iket in set(DetsContainBeta[beta_string]) & set(DetsContainAlpha[alpha_string]):
             *             yield ibra, iket
            */
            // todo: uncomment when all auxiliary arrays are defined
            // m_alpha_singles.foreach_entry_by_key(alpha_channel, [&](uint_t alpha_string){
            //     v_t<FrmOnvField> valid_beta_strings = {};
            //     auto access_beta_with_alpha = m_beta_with_alpha.access(alpha_string);
            //     auto access_beta_singles = m_beta_singles.access(beta_channel);
            //     std::set_intersection(access_beta_with_alpha.m_entry_cbegin, access_beta_with_alpha.m_entry_cend,
            //                           access_beta_singles.m_entry_cbegin, access_beta_singles.m_entry_cend,
            //                           std::back_inserter(valid_beta_strings));
            //     for (auto& beta_string : valid_beta_strings) {
            //         uintv_t dets = {};
            //         auto access_dets_contain_alpha = m_dets_contain_alpha.access(alpha_string);
            //         auto access_dets_contain_beta = m_dets_contain_beta.access(beta_string);
            //         std::set_intersection(access_dets_contain_beta.m_entry_cbegin, access_dets_contain_beta.m_entry_cend,
            //                               access_dets_contain_alpha.m_entry_cbegin, access_dets_contain_alpha.m_entry_cend,
            //                               std::back_inserter(dets));
            //         for (auto& iket : dets){
            //             ket_row.jump(iket);
            //             make_contrib_fn();
            //         }
            //     }
            // });

        }

    void fill_rdm3(PureRdm* rdm) const {
        REQUIRE_TRUE(rdm, "RDM pointer should not be null");
        REQUIRE_TRUE(rdm->m_ranksig == opsig::c_trip, "RDM object should be three-body");
        // todo
    }

    void fill_fock_rdm4(FockRdm4* rdm) const {
        REQUIRE_TRUE(rdm, "RDM pointer should not be null");
        // todo
    }

public:

    /**
     * fill all RDMs
     */
    void fill() {
        {
            auto ptr = m_rdms->get_pure_rdm(opsig::c_sing);
            if (ptr) fill_rdm1(ptr);
        }
        {
            auto ptr = m_rdms->get_pure_rdm(opsig::c_doub);
            if (ptr) fill_rdm2(ptr);
        }
        {
            auto ptr = m_rdms->get_pure_rdm(opsig::c_trip);
            if (ptr) fill_rdm3(ptr);
        }
        {
            auto ptr = m_rdms->m_fock_4rdm;
            if (ptr) fill_fock_rdm4(ptr);
        }
    }


public:
    Caspt2Filler(const Table<MbfWeightRow>& hist, Rdms* rdms):
        m_hist(hist),
        m_fock_x_hist("Fock-perturbed hist WF", MbfWeightRow(hist.m_row),DistribOptions(), Sizing{1000, 1.0}, MbfWeightRow(hist.m_row), Sizing{1000, 1.0}),
        m_rdms(rdms),
        m_work_conns(mbf::get_basis(hist.m_row.m_mbf).size()),
        m_work_com_ops(mbf::get_basis(hist.m_row.m_mbf).size()),
        m_dets_contain_alpha("spin channel to index map (alpha)", hist.m_row.m_mbf.m_format.m_shape[1]),
        m_dets_contain_beta("spin channel to index map (beta)", hist.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_with_alpha("spin channel to spin channel map (alpha)", hist.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_with_beta("spin channel to spin channel map (beta)", hist.m_row.m_mbf.m_format.m_shape[1])
        {

        logging::info("Constructing auxiliary arrays for RDM calculation");
        const auto displ = mpi::evenly_shared_displ(hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(hist.nrow_in_use());
        auto bra = m_hist.m_row;
        buffered::FrmOnvSpinChannel alpha_channel(bra.m_mbf.m_basis.m_nsite);
        buffered::FrmOnvSpinChannel beta_channel(bra.m_mbf.m_basis.m_nsite);

        for (bra.restart(displ); bra.in_range(displ + count); ++bra) {
            bra.m_mbf.copy_alpha_to(alpha_channel);
            bra.m_mbf.copy_beta_to(beta_channel);
            m_dets_contain_alpha.insert(alpha_channel, bra.index());
            m_dets_contain_beta.insert(beta_channel, bra.index());
            m_beta_with_alpha.insert(alpha_channel, beta_channel);
            m_alpha_with_beta.insert(beta_channel, alpha_channel);
        }
        m_dets_contain_alpha.collate();
        m_dets_contain_beta.collate();
        // m_beta_with_alpha.collate();
        // m_alpha_with_beta.collate();

        // {
        //     auto outer_fn = [&](uint_t index, const Smuvi<field::FrmOnvSpinChannel>::SmuviEntriesWithKey &alpha_channel) {
        //         auto tmp_alpha = alpha_channel;
        //         auto inner_fn = [&](uint_t isite) {
        //             tmp_alpha.clr(isite);
        //             m_alpha_single_dict.insert(tmp_alpha, alpha_channel);
        //             tmp_alpha.set(isite);
        //         };
        //         alpha_channel.foreach_setbit(inner_fn);
        //     };
        //     m_dets_contain_alpha.foreach(outer_fn);
        // }

        logging::info("successfully constructed auxiliary arrays for RDM calculation");

        // if there's no Fock*4RDM object allocated, there's nothing left to do
        if (!m_rdms || !m_rdms->m_fock_4rdm) return;
        logging::info("preparing Fock-perturbed vector F |0> from diagonal Fock matrix");
        {
            // if the Fock object is non-diagonal, construct F*psi using the dense matrix overlaod
            auto ptr = dynamic_cast<const NonDiagFockRdm4*>(m_rdms->m_fock_4rdm);
            if (ptr) make_psi1(m_fock_x_hist, ptr->m_fock);
        }
        {
            // if the Fock object is non-diagonal, construct F*psi using the dense vector overlaod
            auto ptr = dynamic_cast<const DiagFockRdm4*>(m_rdms->m_fock_4rdm);
            make_psi1(m_fock_x_hist, ptr->m_fock);
        }
        logging::info("successfully prepared F |0> with {} total rows", mpi::all_sum(m_fock_x_hist.m_store.nrow_in_use()));
    }
};

#endif //M7_CASPT2FILLER_H
