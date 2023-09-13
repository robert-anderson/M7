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

    const conf::Rdms& m_opts;

    communicator::BasicSend<MbfWeightRow, MbfWeightRow> m_psi_pq;
    communicator::BasicSend<MbfWeightRow, MbfWeightRow> m_psi_fock;

    pose::Rdm m_rdm;
    const std::unique_ptr<FockMatrix> m_fock_mat;

    typedef std::pair<std::pair<uint_t, uint_t>, ham_t> pq_val_t;
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

    void make_psi1(communicator::BasicSend<MbfWeightRow, MbfWeightRow>& psi1, const dense::SquareMatrix<ham_t>& mat, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        auto fn = [&pq_vals](uint_t irow, uint_t icol, ham_t elem) {
            pq_vals.emplace_back(uintp_t(irow, icol), elem);
        };
        mat.foreach(fn, tol);
        make_psi1(psi1, pq_vals);
    }

    void fill(const Table<MbfWeightRow>& bra, uint_t bra_displ, uint_t bra_count, const Table<MbfWeightRow>& ket) {
        m_rdm.m_store.clear();
        auto bra_row = bra.m_row;
        auto ket_row = ket.m_row;
        conn::Mbf conn(bra_row.m_mbf);
        for (bra_row.restart(bra_displ); bra_row.in_range(bra_displ + bra_count); ++bra_row) {
            for (ket_row.restart(); ket_row; ++ket_row) {
                auto exsig = mbf::exsig(bra_row.m_mbf, ket_row.m_mbf);
                if (exsig == opsig::c_invalid) continue;
                if (exsig.nfrm_cre() > 2) continue;
                const auto contrib = bra_row.m_weight[0] * ket_row.m_weight[0];
                m_rdm.make_contribs(bra_row.m_mbf, ket_row.m_mbf, contrib);
            }
        }
        m_rdm.communicate();
        auto& recv_row = m_rdm.m_send_recv.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            auto& store_row = m_rdm.m_store.lookup_or_insert(recv_row.m_inds);
            store_row.m_values += recv_row.m_values;
        }
    }

    void fill(const Table<MbfWeightRow>& bra, const Table<MbfWeightRow>& ket) {
        fill(bra, 0ul, bra.nrow_in_use(), ket);
    }

    sys::Sector get_sector() const {
        const auto& row = m_hist.m_row;
        row.restart();
        auto nelec = row.m_mbf.nsetbit();
        return {mbf::get_basis(row.m_mbf), {nelec, {}}};
    }

public:
    Caspt2Filler(const Table<MbfWeightRow>& hist, const conf::Rdms& opts):
        m_hist(hist), m_opts(opts),
        m_psi_pq("excit-perturbed hist WF", MbfWeightRow(hist.m_row),DistribOptions(), Sizing{1000, 1.0}, MbfWeightRow(hist.m_row), Sizing{1000, 1.0}),
        m_psi_fock("Fock-perturbed hist WF", MbfWeightRow(hist.m_row),DistribOptions(), Sizing{1000, 1.0}, MbfWeightRow(hist.m_row), Sizing{1000, 1.0}),
        m_rdm("CASPT2 transition 2RDM", mbf::get_basis(hist.m_row.m_mbf)),
        m_fock_mat(opts.m_fock_4rdm.m_enabled ? new FockMatrix(hist.m_row.m_mbf.m_basis.m_nsite, opts.m_fock_4rdm.m_fock_path) : nullptr){
    }

    void fill_and_save(hdf5::GroupWriter* gw_3300, hdf5::GroupWriter* gw_4400f) {
        const auto nsite = m_hist.m_row.m_mbf.m_basis.m_nsite;
        v_t<pq_val_t> pq_vals(1ul);
        pq_vals.back().second = 1.0;
        for (uint_t p = 0ul; p < nsite; ++p) {
            for (uint_t q = 0ul; q < nsite; ++q) {
                pq_vals.back().first = {p, q};
                logging::info("preparing E_pq |0> with spatial orbital indices p={}, q={}", p, q);
                make_psi1(m_psi_pq, pq_vals);
                logging::info("successfully prepared E_pq |0> with {} total rows", mpi::all_sum(m_psi_pq.m_store.nrow_in_use()));
                if (gw_3300) {
                    fill(m_psi_pq.m_store, m_hist);
                    m_rdm.m_store.save(*gw_3300, logging::format("{},{}", p, q), true);
                }
                if (gw_4400f) {
                    fill(m_psi_pq.m_store, m_psi_fock.m_store);
                    m_rdm.m_store.save(*gw_4400f, logging::format("{},{}", p, q), true);
                }
            }
        }
    }

    void fill_and_save() {

        hdf5::FileWriter fw("M7.pose.h5");
        hdf5::GroupWriter gw_sf(fw, "spinfree");
        hdf5::GroupWriter gw_3300(gw_sf, "3300");

        if (m_fock_mat) {
            hdf5::GroupWriter gw_4400f(gw_sf, "4400f");
            logging::info("preparing F |0> for CASPT2 intermediate");
            make_psi1(m_psi_fock, *m_fock_mat);
            logging::info("successfully prepared F |0> with {} total rows", mpi::all_sum(m_psi_fock.m_store.nrow_in_use()));
            fill_and_save(&gw_3300, &gw_4400f);
        }
        fill_and_save(&gw_3300, nullptr);
    }

};

#endif //M7_CASPT2FILLER_H
