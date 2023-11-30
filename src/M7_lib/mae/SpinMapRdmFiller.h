//
// Created by Robert John Anderson on 22/08/2023.
//

#ifndef M7_SPINMAPRDMFILLER_H
#define M7_SPINMAPRDMFILLER_H

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
class SpinMapRdmFiller {
    /**
     * Histogrammable set of determinants of which to compute the outer product in filling the MAEs
     */
    const Table<MbfWeightRow>& m_bra;
    const Table<MbfWeightRow>& m_ket;

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

    SpinChannelToIndsSmuvi m_dets_contain_alpha;
    SpinChannelToIndsSmuvi m_dets_contain_beta;
    SpinChannelToSpinChannelSmuvi m_beta_with_alpha;
    SpinChannelToSpinChannelSmuvi m_alpha_with_beta;
    SpinChannelToSpinChannelSmuvi m_alpha_single_dict;
    SpinChannelToSpinChannelSmuvi m_beta_single_dict;
    SpinChannelToSpinChannelSmuvi m_alpha_singles;
    SpinChannelToSpinChannelSmuvi m_beta_singles;
    SpinChannelToSpinChannelSmuvi m_alpha_double_dict;
    SpinChannelToSpinChannelSmuvi m_beta_double_dict;
    SpinChannelToSpinChannelSmuvi m_alpha_doubles;
    SpinChannelToSpinChannelSmuvi m_beta_doubles;

    typedef std::pair<std::pair<uint_t, uint_t>, ham_t> pq_val_t;
    /**
     * make a linear combination psi1 = G * psi where G = sum_pq g_pq E_pq and psi is m_hist
     * @param psi1
     *  the output of G * psi
     * @param pq_vals
     *  a sparse representation of the elements of the contraction ((p, q), g_pq)
     */
    static void make_psi1(const Table<MbfWeightRow>& psi0, communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const v_t<pq_val_t>& pq_vals) {
        psi1.m_store.clear();
        const auto displ = mpi::evenly_shared_displ(psi0.nrow_in_use());
        const auto count = mpi::evenly_shared_count(psi0.nrow_in_use());

        auto hist_row = psi0.m_row;
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
    static void make_psi1(const Table<MbfWeightRow>& psi0, communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const dense::SquareMatrix<ham_t>& coeffs, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        auto fn = [&pq_vals](uint_t irow, uint_t icol, ham_t elem) {
            pq_vals.emplace_back(uintp_t(irow, icol), elem);
        };
        coeffs.foreach(fn, tol);
        make_psi1(psi0, psi1, pq_vals);
    }

    /**
     * overload to make psi1 = G * psi where g_ii are elements of a dense vector
     */
    static void make_psi1(const Table<MbfWeightRow>& psi0, communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const dense::Vector<ham_t>& coeffs, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        for (uint_t ielem = 0ul; ielem < coeffs.nelement(); ++ielem)
            pq_vals.emplace_back(uintp_t(ielem, ielem), coeffs[ielem]);
        make_psi1(psi0, psi1, pq_vals);
    }

    // todo: make private again
public:

    void make_contribs(Rdm* rdm, const field::Mbf& src, const field::Mbf& dst, wf_t contrib) const {
        auto& conn = m_work_conns[src];
        auto& com_ops = m_work_com_ops[src];
        conn.connect(src, dst, com_ops);
        rdm->make_contribs(src, conn, com_ops, contrib);
    }

    void fill_rdm(Rdm* rdm) const {
        if (!rdm) return;

        const auto displ = mpi::evenly_shared_displ(m_bra.nrow_in_use());
        const auto count = mpi::evenly_shared_count(m_bra.nrow_in_use());
        auto bra_row = m_bra.m_row;
        auto ket_row = m_ket.m_row;

        auto make_contrib_fn = [&]() {
            const auto contrib = bra_row.m_weight[0] * ket_row.m_weight[0];
            make_contribs(rdm, bra_row.m_mbf, ket_row.m_mbf, contrib);
        };

        const auto order_fn = [&](const field::FrmOnvSpinChannel& i, const field::FrmOnvSpinChannel& j) -> bool {
            return i < j;
        };

        buffered::FrmOnvSpinChannel alpha_channel(bra_row.m_mbf.m_basis.m_nsite);
        buffered::FrmOnvSpinChannel beta_channel(bra_row.m_mbf.m_basis.m_nsite);
        uint_t counter = 0;
        for (bra_row.restart(displ); bra_row.in_range(displ + count); ++bra_row) {
            counter += 1;
            if (counter%1000 == 0) logging::info("currently in iteration {}", counter);
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

            if (rdm->m_ranksig == opsig::c_sing) continue;

            // alpha-beta
            m_alpha_singles.foreach_value(alpha_channel, [&](const field::FrmOnvSpinChannel &alpha_string){
                const auto indices_dets_with_alpha = m_dets_contain_alpha.access(alpha_string);
                m_beta_with_alpha.foreach_common_value(alpha_string, m_beta_singles, beta_channel,
                                                      [&](const field::FrmOnvSpinChannel &common_string){
                    indices_dets_with_alpha.m_value_row.jump(common_string.m_row->index());
                    const uint_t iket = indices_dets_with_alpha.m_value_row.m_value;
                    ket_row.jump(iket);
                    DEBUG_ASSERT_EQ(bra_row.m_mbf.nbeta_not_in(ket_row.m_mbf), 1, "only beta singles yield valid contributions.");
                    make_contrib_fn();
                }, order_fn);
            });

            if (rdm->m_ranksig == opsig::c_doub) continue;

            // 4x(alpha) 2x(beta), here alpha_doubles instead of alpha_singles, otherwise exact copy of 2RDM code
            m_alpha_doubles.foreach_value(alpha_channel, [&](const field::FrmOnvSpinChannel &alpha_string){
               const auto indices_dets_with_alpha = m_dets_contain_alpha.access(alpha_string);
               m_beta_with_alpha.foreach_common_value(alpha_string, m_beta_singles, beta_channel,
                                                     [&](const field::FrmOnvSpinChannel &common_string){
                   indices_dets_with_alpha.m_value_row.jump(common_string.m_row->index());
                   const uint_t iket = indices_dets_with_alpha.m_value_row.m_value;
                   ket_row.jump(iket);
                   DEBUG_ASSERT_EQ(bra_row.m_mbf.nbeta_not_in(ket_row.m_mbf), 1, "only beta singles yield valid contributions.");
                   make_contrib_fn();
               }, order_fn);
            });
            // 2x(alpha) 4x(beta), flip the roles of alpha and beta
            m_beta_doubles.foreach_value(beta_channel, [&](const field::FrmOnvSpinChannel &beta_string){
               const auto indices_dets_with_beta = m_dets_contain_beta.access(beta_string);
               m_alpha_with_beta.foreach_common_value(beta_string, m_alpha_singles, alpha_channel,
                                                     [&](const field::FrmOnvSpinChannel &common_string){
                   indices_dets_with_beta.m_value_row.jump(common_string.m_row->index());
                   const uint_t iket = indices_dets_with_beta.m_value_row.m_value;
                   ket_row.jump(iket);
                   DEBUG_ASSERT_EQ(bra_row.m_mbf.nalpha_not_in(ket_row.m_mbf), 1, "only alpha singles yield valid contributions.");
                   make_contrib_fn();
               }, order_fn);
            });
        }
    }


public:
    SpinMapRdmFiller(const Table<MbfWeightRow>& hist_bra, const Table<MbfWeightRow>& hist_ket):
        m_bra(hist_bra),
        m_ket(hist_ket),
        m_work_conns(mbf::get_basis(m_ket.m_row.m_mbf).size()),
        m_work_com_ops(mbf::get_basis(m_ket.m_row.m_mbf).size()),
        m_dets_contain_alpha("spin channel to index map (alpha)", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_dets_contain_beta("spin channel to index map (beta)", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_with_alpha("spin channel to spin channel map (alpha)", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_with_beta("spin channel to spin channel map (beta)", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_single_dict("auxiliary spin channel to spin channel map alpha singles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_single_dict("auxiliary spin channel to spin channel map beta singles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_singles("spin channel to spin channel map alpha singles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_singles("spin channel to spin channel map beta singles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_double_dict("auxiliary spin channel to spin channel map alpha doubles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_double_dict("auxiliary spin channel to spin channel map beta doubles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_alpha_doubles("spin channel to spin channel map alpha doubles", m_ket.m_row.m_mbf.m_format.m_shape[1]),
        m_beta_doubles("spin channel to spin channel map beta doubles", m_ket.m_row.m_mbf.m_format.m_shape[1]) {

        logging::info("Constructing auxiliary arrays for RDM calculation");
        const auto displ = mpi::evenly_shared_displ(hist_ket.nrow_in_use());
        const auto count = mpi::evenly_shared_count(hist_ket.nrow_in_use());
        auto hist_row = m_ket.m_row;

        buffered::FrmOnvSpinChannel alpha_channel(hist_row.m_mbf.m_basis.m_nsite);
        buffered::FrmOnvSpinChannel beta_channel(hist_row.m_mbf.m_basis.m_nsite);
        /**
         *  Construct SMUVIs which given a FrmOnvSpinChannel yield
         *      the rows of the histogrammed set containing this spin string: m_dets_contain_(spin),
         *      the opposite spin strings occurring with it: m_(spin1)_with_(spin2)
         */
        for (hist_row.restart(displ); hist_row.in_range(displ + count); ++hist_row) {
            hist_row.m_mbf.copy_alpha_to(alpha_channel);
            hist_row.m_mbf.copy_beta_to(beta_channel);
            m_dets_contain_alpha.insert(alpha_channel, hist_row.index());
            m_dets_contain_beta.insert(beta_channel, hist_row.index());
            m_beta_with_alpha.insert(alpha_channel, beta_channel);
            m_alpha_with_beta.insert(beta_channel, alpha_channel);
        }
        /**
         * m_dets_contain_alpha/beta are (bit string -> det index) maps, but the indices need to be ordered by the
         * complementary spin string (beta for alpha) of the determinants which are encoded by the indices.
         */
        {
            // these serve merely as temporary buffers for the comparison function
            buffered::FrmOnvSpinChannel beta_channel1(hist_row.m_mbf.m_basis.m_nsite);
            buffered::FrmOnvSpinChannel beta_channel2(hist_row.m_mbf.m_basis.m_nsite);
            const auto order_fn = [&](const field::Number<uint_t>& i, const field::Number<uint_t>& j) -> bool {
                m_ket.m_row.jump(i);
                m_ket.m_row.m_mbf.copy_beta_to(beta_channel1);
                m_ket.m_row.jump(j);
                m_ket.m_row.m_mbf.copy_beta_to(beta_channel2);
                return beta_channel1 < beta_channel2;
            };
            m_dets_contain_alpha.collate(order_fn);
        }
        {
            // same as for beta
            buffered::FrmOnvSpinChannel alpha_channel1(hist_row.m_mbf.m_basis.m_nsite);
            buffered::FrmOnvSpinChannel alpha_channel2(hist_row.m_mbf.m_basis.m_nsite);
            const auto order_fn = [&](const field::Number<uint_t>& i, const field::Number<uint_t>& j) -> bool {
                m_ket.m_row.jump(i);
                m_ket.m_row.m_mbf.copy_alpha_to(alpha_channel1);
                m_ket.m_row.jump(j);
                m_ket.m_row.m_mbf.copy_alpha_to(alpha_channel2);
                return alpha_channel1 < alpha_channel2;
            };
            m_dets_contain_beta.collate(order_fn);
        }
        // unlike for (bit string -> det index), (bit string -> bit string) maps can reuse the same order_fn
        const auto order_fn = [&](const field::FrmOnvSpinChannel& i, const field::FrmOnvSpinChannel& j) -> bool {
            return i < j;
        };
        m_alpha_with_beta.collate(order_fn);
        m_beta_with_alpha.collate(order_fn);
        /**
         *  Generate all (N - 1) electron states from the spin strings.
         */
        SpinChannelToSpinChannelSmuvi* spin_single_dict = nullptr;
        buffered::FrmOnvSpinChannel tmp_spin_channel(hist_row.m_mbf.m_basis.m_nsite);
        auto gen_one_less_electron= [&](const field::FrmOnvSpinChannel& key, SpinChannelToIndsSmuvi::AccessResult idets){
            tmp_spin_channel = key;
            auto inner_fn = [&](uint_t isite) {
                tmp_spin_channel.clr(isite);
                spin_single_dict->insert(tmp_spin_channel, key);
                tmp_spin_channel.set(isite);
            };
            key.foreach_setbit(inner_fn);
        };
        spin_single_dict = &m_alpha_single_dict;
        m_dets_contain_alpha.foreach_key(gen_one_less_electron);
        m_alpha_single_dict.collate(order_fn);
        spin_single_dict = &m_beta_single_dict;
        m_dets_contain_beta.foreach_key(gen_one_less_electron);
        m_beta_single_dict.collate(order_fn);
        /**
         *  Loop over all (N - 1) electron keys of the m_(spin)_single_dict SMUVI and add pairs as key value pairs into
         *  the m_(spin)_singles SMUVIs. For example, the key [0011100] may point to [1011100, 0111100, 0011101, ...], then
         *      m_(spin)_singles[1011100] = [0111100, 0011101, ...]
         *      m_(spin)_singles[0111100] = [1011100, 0011101, ...]
         *      m_(spin)_singles[0011101] = [1011100, 0111100, ...]
         *   All these determinants are single excitations from each other.
         */
        SpinChannelToSpinChannelSmuvi* spin_singles = nullptr;
        auto gen_spin_singles = [&](const field::FrmOnvSpinChannel& key, SpinChannelToSpinChannelSmuvi::AccessResult strings){
             spin_single_dict->foreach_value_pair(key, [&](const field::FrmOnvSpinChannel& value_1, const field::FrmOnvSpinChannel& value_2){
                 spin_singles->insert(value_1, value_2);
             });
        };
        spin_single_dict = &m_alpha_single_dict;
        spin_singles = &m_alpha_singles;
        m_alpha_single_dict.foreach_key(gen_spin_singles);
        m_alpha_singles.collate(order_fn);
        spin_single_dict = &m_beta_single_dict;
        spin_singles = &m_beta_singles;
        m_beta_single_dict.foreach_key(gen_spin_singles);
        m_beta_singles.collate(order_fn);
        /**
         *   Generate all (N - 2) electron states from the spin strings.
         */
        SpinChannelToSpinChannelSmuvi* spin_double_dict = nullptr;
        auto gen_two_less_electron= [&](const field::FrmOnvSpinChannel& key, SpinChannelToIndsSmuvi::AccessResult idets){
            tmp_spin_channel = key;
            auto inner_fn = [&](uint_t isite1, uint_t isite2) {
                tmp_spin_channel.clr(isite1);
                tmp_spin_channel.clr(isite2);
                spin_double_dict->insert(tmp_spin_channel, key);
                tmp_spin_channel.set(isite1);
                tmp_spin_channel.set(isite2);
            };
            key.foreach_setbit_pair(inner_fn);
        };
        spin_double_dict = &m_alpha_double_dict;
        m_dets_contain_alpha.foreach_key(gen_two_less_electron);
        m_alpha_double_dict.collate(order_fn);
        spin_double_dict = &m_beta_double_dict;
        m_dets_contain_beta.foreach_key(gen_two_less_electron);
        m_beta_double_dict.collate(order_fn);
        /**
         * Loop over each key of the (N - 2) electron SMUVI and add the values to the m_(spin)_doubles SMUVI analogous
         * to the singles, but take care that only genuine doubles are counted.
         */
        SpinChannelToSpinChannelSmuvi* spin_doubles = nullptr;
        auto gen_spin_doubles = [&](const field::FrmOnvSpinChannel& key, SpinChannelToSpinChannelSmuvi::AccessResult strings){
            spin_double_dict->foreach_value_pair(key, [&](const field::FrmOnvSpinChannel& value_1, const field::FrmOnvSpinChannel& value_2){
                if (value_1.nsetbit_not_in(value_2) == 2) spin_doubles->insert(value_1, value_2);
            });
        };
        spin_double_dict = &m_alpha_double_dict;
        spin_doubles = &m_alpha_doubles;
        m_alpha_double_dict.foreach_key(gen_spin_doubles);
        m_alpha_doubles.collate(order_fn);
        spin_double_dict = &m_beta_double_dict;
        spin_doubles = &m_beta_doubles;
        m_beta_double_dict.foreach_key(gen_spin_doubles);
        m_beta_doubles.collate(order_fn);

        m_dets_contain_alpha.remap_accessors();
        m_dets_contain_beta.remap_accessors();
        m_beta_with_alpha.remap_accessors();
        m_alpha_with_beta.remap_accessors();
        m_alpha_single_dict.remap_accessors();
        m_beta_single_dict.remap_accessors();
        m_alpha_singles.remap_accessors();
        m_beta_singles.remap_accessors();
        m_alpha_double_dict.remap_accessors();
        m_beta_double_dict.remap_accessors();
        m_alpha_doubles.remap_accessors();
        m_beta_doubles.remap_accessors();

        logging::info("successfully constructed auxiliary arrays for RDM calculation");
    }


    /**
     * fill all RDMs
     */
    static void fill(const Table<MbfWeightRow>& hist, Rdms* rdms) {
        if (!rdms) return;

        // TODO: how to avoid code duplication?
        auto& row = hist.m_row;
        wf_comp_t norm = 0.0;
        for (row.restart(); row; ++row) norm += math::pow<2>(std::abs(row.m_weight[0]));
        if (mpi::i_am_root()) rdms->m_total_norm.m_local = norm;

        bool have_pure = false;
        have_pure |= rdms->get_pure_rdm(opsig::c_sing) != nullptr;
        have_pure |= rdms->get_pure_rdm(opsig::c_doub) != nullptr;
        have_pure |= rdms->get_pure_rdm(opsig::c_trip) != nullptr;
        if (have_pure) {
            // at least one of the pure RDM instances is allocated, so make aux arrays for the hist-hist RDMs and fill
            SpinMapRdmFiller filler(hist, hist);
            {
                auto ptr = rdms->get_pure_rdm(opsig::c_sing);
                if (ptr) filler.fill_rdm(ptr);
            }
            {
                auto ptr = rdms->get_pure_rdm(opsig::c_doub);
                if (ptr) filler.fill_rdm(ptr);
            }
            {
                auto ptr = rdms->get_pure_rdm(opsig::c_trip);
                if (ptr) filler.fill_rdm(ptr);
            }
        }
        if (rdms->m_fock_4rdm) {
            /*
             * the result of F * m_hist where F = sum_pq f_pq E_pq
             */
            communicator::BasicSend<MbfWeightRow, MbfWeightRow> fock_x_hist(
                    "Fock-perturbed hist WF", MbfWeightRow(hist.m_row),DistribOptions(), Sizing{1000, 1.0}, MbfWeightRow(hist.m_row), Sizing{1000, 1.0});

            logging::info("preparing Fock-perturbed vector F |0> from diagonal Fock matrix");
            {
                // if the Fock object is non-diagonal, construct F*psi using the matrix overload
                auto ptr = dynamic_cast<const NonDiagFockRdm4*>(rdms->m_fock_4rdm);
                if (ptr) make_psi1(hist, fock_x_hist, ptr->m_fock);
            }
            {
                // if the Fock object is diagonal, construct F*psi using the vector overload
                auto ptr = dynamic_cast<const DiagFockRdm4 *>(rdms->m_fock_4rdm);
                if (ptr) make_psi1(hist, fock_x_hist, ptr->m_fock);
            }
            logging::info("successfully prepared F |0> with {} total rows", mpi::all_sum(fock_x_hist.m_store.nrow_in_use()));

            SpinMapRdmFiller(hist, fock_x_hist.m_store).fill_rdm(rdms->m_fock_4rdm);
        }
    }
};

#endif //M7_SPINMAPRDMFILLER_H
