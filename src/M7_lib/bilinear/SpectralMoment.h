//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_SPECTRALMOMENT_H
#define M7_SPECTRALMOMENT_H

#include <M7_lib/mae/MaeTable.h>
#include <M7_lib/conf/Conf.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/util/Vector.h>
#include "M7_lib/propagator/ExactLinear.h"



//Communicator<MaeRow<wf_t>, MaeRow<wf_t>, true>
//struct SpecMom {
//    const uint_t m_order;
//    SendRecv<>
//    SpecMom(uint_t order): m_exsig(exsig), m_order(order){
//        REQUIRE_EQ(order, 1ul, "Spectral moment quantities are currently only implemented for n=1")
//        REQUIRE_TRUE(m_exsig.is_pure_frm(), "Spectral moment excitations must refer to purely fermionic perturbations");
//    }
//};

/**
 * the classes defined here relate to the estimation of the following kinds of expectation values:
 * order-n spectral moment of the hole kind:
 *  <Psi | p+ H^n q | Psi>
 * order-n spectral moment of the particle kind:
 *  <Psi | p H^n q+ | Psi>
 */
namespace spec_mom {
    /**
     * Row type is generic to hole/particle kinds
     */
    struct CommRow : public Row {
        field::Mbf m_dst_mbf;
        field::Number<wf_t> m_hnci_factor;
        field::Flag m_src_deterministic;
        field::Number<uint8_t> m_ipart_dst;
        field::Number<uint8_t> m_icre;
        field::Number<uint8_t> m_iann;

        explicit CommRow(const sys::Basis& basis) :
            m_dst_mbf(this, basis, "destination MBF"),
            m_hnci_factor(this, "source weight multiplied by H matrix elements and with the correct probabilistic normalization"),
            m_src_deterministic(this, "source deterministic flag"),
            m_ipart_dst(this, "WF part index of destination"),
            m_icre(this, "index of the spin orbital created in the expectation value"),
            m_iann(this, "index of the spin orbital annihilated in the expectation value"){}
    };

    typedef SendRecv<CommRow, Table<CommRow>> comm_t;

    struct SpecMom {
        const uint_t m_order;
        const bool m_hole_kind;
        /**
         * communication tables
         */
        comm_t m_comm;
        /**
         * the sampled quantities only have a memory requirement of creation values
         */
        dense::Matrix<ham_t> m_matrix;

        SpecMom(uint_t order, bool hole_kind, const Wavefunction& wf, uintv_t cres, uintv_t anns):
                m_order(order), m_hole_kind(hole_kind), m_wf(wf),
                m_frm_cres(std::move(cres)), m_frm_anns(std::move(std::move(anns))),
                m_icre_to_irow(make_to_row_or_col(m_frm_cres, m_wf.m_sector.basis().m_frm.m_nspinorb)),
                m_iann_to_icol(make_to_row_or_col(m_frm_anns, m_wf.m_sector.basis().m_frm.m_nspinorb)),
                m_comm(logging::format("order {} {} spec mom", m_order, m_hole_kind ? "hole": "particle"),
                   CommRow(m_wf.m_sector.basis()), {1000, 2}),
                m_matrix(m_frm_cres.size(), m_frm_anns.size()){}

        /**
         * assume all spin orbitals are included
         */
        SpecMom(uint_t order, bool hole_kind, const Wavefunction& wf):
                SpecMom(order, hole_kind, wf,
                        vector::range(0ul, m_wf.m_sector.basis().m_frm.m_nspinorb, 1ul),
                        vector::range(0ul, m_wf.m_sector.basis().m_frm.m_nspinorb, 1ul)){}

        CommRow& add_send_contrib(
            const Mbf& dst_mbf,
            uint_t icre,
            uint_t iann,
            ham_t hnci_factor,
            bool src_deterministic,
            uint_t ipart_dst)
        {
            // if this is a hole moment, the creation operator must be occupied in dst MBF
            DEBUG_ASSERT_TRUE(!m_hole_kind || dst_mbf.get(icre), "this contrib is zero since creation op destroys the bra");
            // if this is a particle moment, the annihilation operator must be vacant in dst MBF
            DEBUG_ASSERT_TRUE(m_hole_kind || !dst_mbf.get(icre), "this contrib is zero since annihilation op destroys the bra");
            // the sampled element must be one specified by the user:
            DEBUG_ASSERT_NE(m_icre_to_irow[icre], ~0ul, "not a sampled creation operator index");
            DEBUG_ASSERT_NE(m_iann_to_icol[iann], ~0ul, "not a sampled annihilation operator index");

            auto& dst_table = m_comm.send(m_wf.m_dist.irank(dst_mbf));

            auto& send_row = dst_table.m_row;
            send_row.push_back_jump();

            send_row.m_dst_mbf = dst_mbf;
            send_row.m_hnci_factor = hnci_factor;
            send_row.m_src_deterministic = src_deterministic;
            send_row.m_ipart_dst = ipart_dst;
            send_row.m_icre = icre;
            send_row.m_iann = iann;
            return send_row;
        }
    };

    struct SpecMoms {
        const uint_t m_max_order;
        const Wavefunction& m_wf;
        /**
         * allow the user to only probe a block of the spectral moment matrix
         */
        const uintv_t m_frm_cres;
        const uintv_t m_frm_anns;
        /**
         * maps from the spinorbs to the row/col indices of the matrix
         */
        const uintv_t m_icre_to_irow;
        const uintv_t m_iann_to_icol;

        conn::Mbf m_conn_work;

        v_t<SpecMom> m_hole;
        v_t<SpecMom> m_particle;

    private:
        static uintv_t make_to_row_or_col(const uintv_t& inds, uint_t nspinorb) {
            // initially, all are invalid
            uintv_t map(nspinorb, ~0ul);
            for (uint_t i=0ul; i<inds.size(); ++i) map[inds[i]] = i;
        }

        /**
         * loop over all given spin-orbital indices and remove from mbf before calling the given functor
         * of course mbf will be modified each iteration, but we can promise constness since the mbf is as returned to
         * the initial state upon return
         */
        template<typename fn_t>
        static void foreach_occ(const Mbf& mbf, const uintv_t ispinorbs, const fn_t& fn) {
            functor::assert_prototype<void(uint_t /*ispinorb*/)>(fn);
            auto& nc_mbf = const_cast<Mbf&>(mbf);
            for (auto ispinorb: ispinorbs) {
                // skip when the spin orb is already vacant
                if (!mbf.get(ispinorb)) continue;
                // now set the mbf to the required state
                nc_mbf.clr(ispinorb);
                // call an arbitrary function
                fn(ispinorb);
                // reset the mbf to its state at input
                nc_mbf.set(ispinorb);
            }
        }

        /**
         * as above but for inserting particles rather than holes
         */
        template<typename fn_t>
        static void foreach_vac(const Mbf& mbf, const uintv_t ispinorbs, const fn_t& fn) {
            functor::assert_prototype<void(uint_t /*ispinorb*/)>(fn);
            auto& nc_mbf = const_cast<Mbf&>(mbf);
            for (auto ispinorb: ispinorbs) {
                // skip when the spin orb is already occupied
                if (!mbf.get(ispinorb)) continue;
                // now set the mbf to the required state
                nc_mbf.set(ispinorb);
                // call an arbitrary function
                fn(ispinorb);
                // reset the mbf to its state at input
                nc_mbf.clr(ispinorb);
            }
        }

    public:
        //SpecMoms(const conf::SpecMoms& /*opts*/) {}
        SpecMoms(uint_t max_order, const Wavefunction& wf, uintv_t cres, uintv_t anns):
            m_max_order(max_order), m_wf(wf),
            m_frm_cres(std::move(cres)), m_frm_anns(std::move(std::move(anns))),
        m_icre_to_irow(make_to_row_or_col(m_frm_cres, m_wf.m_sector.basis().m_frm.m_nspinorb)),
        m_iann_to_icol(make_to_row_or_col(m_frm_anns, m_wf.m_sector.basis().m_frm.m_nspinorb)),
        m_src_work(wf.m_sector), m_conn_work(m_src_work){}
    private:
        void add_all_to_send_exact(const ExactLinear& prop, const Mbf& src, uint_t n, uint_t ispinorb, bool hole, ham_t hnci_factor, bool any_offdiag) {
            if (n > m_max_order) return;
            // working MBF for destinations
            buffered::Mbf dst(src);
            //m_ham.m_frm.make_foreach_iters()
            auto fn = [&]() -> void {
                m_conn_work.apply(src, dst);
                const auto helem = prop.m_ham.get_element(m_src_work, m_conn_work);
                /*
                 * add to current order (n) matrices
                 */
                if (hole) {
                    for (auto icre: m_frm_cres){
                        // skip creation indices which would destroy the dst state
                        if (!dst.get(icre)) continue;
                        m_hole[n].add_send_contrib(dst, icre, ispinorb, hnci_factor*helem, false, 0ul);
                    }
                }
                if (hole) {
                    for (auto icre: m_frm_cres){
                        m_hole[n].add_send_contrib(dst, icre, ispinorb, hnci_factor*helem, false, 0ul);
                    }
                }
                add_all_to_send_exact(prop, dst, n+1, hnci_factor*helem, any_offdiag);
            };
            prop.conn_iters().loop(m_conn_work, m_src_work, fn);
        }

    public:

        /**
         * loop over all holes and particles that can validly be inserted into the
         */
        void add_all_to_send(const ExactLinear& prop, const Mbf& src, wf_t weight) {
            foreach_occ(src, m_frm_anns,
                        [&](uint_t iann) { add_all_to_send_exact(prop, src, 0, iann, true, weight, false); });
            foreach_vac(src, m_frm_cres,
                        [&](uint_t icre) { add_all_to_send_exact(prop, src, 0, icre, false, weight, false); });
        }

//        void add_all_to_send(const Mbf& src, const ExactLinear& prop) {
//            for (auto icre: m_frm_cres)
//
//        }
    };
}

class SpecMoms {

public:


    operator bool() const {
        return false;//!m_active_ranksigs.empty();
    }

    bool all_stores_empty() const {
//        for (auto& ranksig: m_active_ranksigs)
//            if (!m_rdms[ranksig]->m_store.is_freed())
//                return false;
            return true;
    }
};

#endif //M7_SPECTRALMOMENT_H
