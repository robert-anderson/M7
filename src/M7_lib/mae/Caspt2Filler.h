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

#if 0
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

    communicator::BasicSend<MbfWeightRow, MbfWeightRow> m_psi1;

    PureRdm m_rdm;
    SpinFreeRdm m_sf_rdm;

    typedef std::pair<std::pair<uint_t, uint_t>, ham_t> pq_val_t;
    void make_psi1(const v_t<pq_val_t>& pq_vals) {
        m_psi1.m_store.clear();
        const auto displ = mpi::evenly_shared_displ(m_hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(m_hist.nrow_in_use());

        auto hist_row = m_hist.m_row;
        buffered::Mbf dst(hist_row.m_mbf);
        v_t<std::pair<conn::Mbf, ham_t>> conn_vals;
        for (auto& pq_val: pq_vals) {
            const auto p = pq_val.first.first;
            const auto q = pq_val.first.second;
            const auto contract_val = pq_val.second;
            for (uint_t ispin = 0ul; ispin < 2ul; ++ispin) {
                conn_vals.emplace_back(dst, contract_val);
                conn_vals.back().first.m_cre.set(dst.m_basis.ispinorb(ispin, p));
                conn_vals.back().first.m_ann.set(dst.m_basis.ispinorb(ispin, q));
            }
        }

        for (hist_row.restart(displ); hist_row.in_range(displ + count); ++hist_row) {
            for (auto& conn_val: conn_vals) {
                const auto& conn = conn_val.first;
                if (mbf::destroys(conn, dst)) continue;
                const auto& val = conn_val.second;
                conn.apply(hist_row.m_mbf, dst);
                const auto phase = conn.phase(hist_row.m_mbf);
                auto irank_dst = m_psi1.m_dist.irank(dst);
                auto& send_row = m_psi1.m_send_recv.send(irank_dst).m_row;
                send_row.push_back_jump();
                send_row.m_mbf = dst;
                send_row.m_weight = hist_row.m_weight;
                send_row.m_weight *= conn_val.second;
                if (phase) send_row.m_weight *= -1.0;
            }
        }

        m_psi1.communicate();
        /*
         * do mini annihilation loop
         */
        auto& recv_row = m_psi1.m_send_recv.recv().m_row;
        for (recv_row.restart(); recv_row; ++recv_row) {
            auto& dst = m_psi1.m_store.lookup_or_insert(recv_row.m_mbf);
            dst.m_weight += recv_row.m_weight;
        }
    }

    void make_psi1(const dense::SquareMatrix<ham_t>& mat, ham_comp_t tol=1e-10) {
        v_t<pq_val_t> pq_vals;
        auto fn = [&pq_vals](uint_t irow, uint_t icol, ham_t elem) {
            pq_vals.emplace_back(uintp_t(irow, icol), elem);
        };
        mat.foreach(fn, tol);
        make_psi1(pq_vals);
    }

    void fill(const Table<MbfWeightRow>& bra, uint_t bra_displ, uint_t bra_count, const Table<MbfWeightRow>& ket) {
        m_rdm.m_store.clear();
        auto bra_row = bra.m_row;
        auto ket_row = ket.m_row;
        conn::Mbf conn(bra_row.m_mbf);
        buffered::RdmInds rdm_inds(m_rdm.m_ranksig);
        for (bra_row.restart(bra_displ); bra_row.in_range(bra_displ + bra_count); ++bra_row) {
            for (ket_row.restart(); ket_row; ++ket_row) {
                auto exsig = mbf::exsig(bra_row.m_mbf, ket_row.m_mbf);
                if (exsig == opsig::c_invalid) continue;
                if (exsig.nfrm_cre() > 2) continue;
                conn.connect(bra_row.m_mbf, ket_row.m_mbf);
                const auto phase = conn.phase(bra_row.m_mbf);
                const auto contrib = bra_row.m_weight[0] * ket_row.m_weight[0];
                m_rdm.make_full_contrib(rdm_inds, opsig::c_invalid, contrib, phase);
            }
        }
    }

public:
    Caspt2Filler(const Table<MbfWeightRow>& hist, const conf::Rdms& opts, const sys::Particles& particles):
        m_hist(hist), m_psi1("perturbed hist WF", MbfWeightRow(hist.m_row), DistribOptions(), Sizing{}, MbfWeightRow(hist.m_row), Sizing{}),
        m_rdm(opts, opsig::c_doub, sys::Sector(mbf::get_basis(hist.m_row.m_mbf), particles), 1, "CASPT2 transition 2RDM"),
        m_sf_rdm(){
    }

    void fill_and_save() {

    }

};
#endif //M7_CASPT2FILLER_H

#endif //M7_CASPT2FILLER_H
