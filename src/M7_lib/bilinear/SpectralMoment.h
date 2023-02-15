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
 *  An(-1)[p, q] = <Psi | p+ H^n q | Psi>
 * order-n spectral moment of the particle kind:
 *  An(+1)[p, q] = <Psi | p H^n q+ | Psi>
 *
 * contributions to An(+/-1)[p, q] come from a source (bra) N-electron MBF "src", a destination (ket) N-electron MBF "dst"
 * via a chain of n (N+/-1)-electron MBFs (k0, k1, k2, ..., k_n)
 *
 * these are categorised according to the following scheme:
 *  - DIAGONAL: p == q, dst == src, ki == q |src> if hole kind, q+ |src> if particle kind => accounted for by walker
 *      death and block averaging
 *  - OFF-DIAGONAL: all other contributions. note that contribution diagonality refers to the walkers, Hamiltonian
 *      matrix elements, *and* moment indices involved, so even when p == q, dst == src but the k-chain involves
 *      off-diagonal elements, the contributions are considered to be off-diagonal
 *      - DETERMINISTIC-OFF-DIAGONAL:
 *         when src and dst are in the deterministic subspace, and all ki are in the deterministic hole / particle bases
 *      - NON-DETERMINISTIC-OFF-DIAGONAL:
 *         all other non-diagonal contributions
 *
 * the shorthand used for these categories is Diag, DetermOffDiag, NonDetermOffDiag
 */

/**
 * the spectral moment analogue of the Propagator class - responsible for generating contributions of the
 * non-deterministic off-diagonal kind, discarding all those values of (p, q, src, dst, ki) which do not conform to
 * this kind of contribution
 */



struct SpecMom {
    /**
     * spin-orbital indices indexing the perturbers (a subset of all spinorbs)
     */
    const uintv_t m_selected_spinorbs;
    /**
     * the sampled quantities only have a memory requirement of creation values
     */
    v_t<ham_t> m_store;
    /**
     * nspinorb element vector mapping perturber spinorbs gaplessly on [0, number of selected spinorbs)
     */
    const uintv_t m_spinorb_to_store_ind;

private:
    static uintv_t make_spinorb_to_store_ind(const uintv_t& selected_spinorbs, uint_t nspinorb) {
        uintv_t tmp(nspinorb, ~0ul);
        uint_t i = 0ul;
        for (auto ispinorb : selected_spinorbs) tmp[ispinorb] = i++;
        return tmp;
    }

    uint_t store_ind(uint_t ispinorb_left, uint_t ispinorb_right) const {
        auto ileft = m_spinorb_to_store_ind[ispinorb_left];
        DEBUG_ASSERT_NE(ileft, ~0ul, "this spinorb is not a selected perturber");
        auto iright = m_spinorb_to_store_ind[ispinorb_right];
        DEBUG_ASSERT_NE(iright, ~0ul, "this spinorb is not a selected perturber");
        return integer::trigmap_unordered(ileft, iright);
    }

public:

    SpecMom(uintv_t selected_spinorbs, uint_t nspinorb):
        m_selected_spinorbs(std::move(selected_spinorbs)),
        m_spinorb_to_store_ind(make_spinorb_to_store_ind(m_selected_spinorbs, nspinorb)){}

    void average() {
        auto recv = m_store;
        mpi::all_sum(m_store.data(), recv.data(), recv.size());
        m_store = recv;
    }

    void make_contrib(uint_t ispinorb_left, uint_t ispinorb_right, ham_t contrib) {
        m_store[store_ind(ispinorb_left, ispinorb_right)] += contrib;
    }
};

//typedef SendRecv<DstFinderRow, Table<DstFinderRow>> dst_finder_comm_t;
///**
// * communication tables
// */
//dst_finder_comm_t m_comm;
#if 0

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

#endif //M7_SPECTRALMOMENT_H


struct SpecMomGen {

};
struct ExactSpecMomGen {

};
struct StochExactSpecMomGen {

};



struct SpecMoms {
    /**
     * Row type is generic to hole/particle kinds. When a connection is generated, we must communicate the contribution
     * to the MPI rank that stores the instantaneous value of the destination MBF
     */
    struct DstFinderRow : public Row {
        field::Mbf m_dst_mbf;
        field::Number<wf_t> m_hnci_factor;
        field::Number<uint8_t> m_ham_order;
        field::Number<uint8_t> m_ipart_dst;
        field::Number<uint8_t> m_ispinorb_left;
        field::Number<uint8_t> m_ispinorb_right;

        explicit DstFinderRow(const sys::Basis& basis) :
            m_dst_mbf(this, basis, "destination MBF"),
            m_hnci_factor(this, "source weight multiplied by H matrix elements and with the correct probabilistic normalization"),
            m_ham_order(this, "power of the Hamiltonian in the estimated expectation value"),
            m_ipart_dst(this, "WF part index of destination"),
            m_ispinorb_left(this, "index of the perturbing spin orbital at the left of the moment expectation value"),
            m_ispinorb_right(this, "index of the perturbing spin orbital at the right of the moment expectation value"){}
    };
    typedef SendRecv<DstFinderRow, Table<DstFinderRow>> dst_finder_comm_t;
    /**
     * const ref to wavefunction class so we can lookup dst walker weights
     */
    const Wavefunction& m_wf;
    /**
     * ref to propagator class which is used to generate connections in the perturbed many-body basis
     */
    Propagator& m_prop;
    /**
     * maximum power of the Hamiltonian for which to compute the moments
     */
    const uint_t m_max_order;
    /**
     * communication tables for finding the dst weights in hole and particle-type moments
     */
    dst_finder_comm_t m_hole_dst_finder_table;
    dst_finder_comm_t m_particle_dst_finder_table;
    /**
     * user-definable subset of spinorbs to comprise the set of perturbers
     */
    const uintv_t m_selected_spinorbs;
    /**
     * hole and particle vectors of average moment-accumulating objects, one for each power of H
     */
    v_t<SpecMom> m_hole_spec_moms;
    v_t<SpecMom> m_particle_spec_moms;


private:
    static uintv_t make_selected_spinorbs(const conf::SpecMoms& opts, uint_t nspinorb){
        if (opts.m_spinorbs.m_value.empty()) return vector::range(0ul, nspinorb);
        return opts.m_spinorbs;
    }

    DstFinderRow& add_send_contrib(
        bool hole,
        const Mbf& dst_mbf,
        ham_t hnci_factor,
        uint_t ham_order,
        uint_t ipart_dst,
        uint_t ispinorb_left,
        uint_t ispinorb_right)
    {
        // if this is a hole moment contribution, the creation operator must be occupied in dst MBF
        DEBUG_ASSERT_TRUE(!hole || dst_mbf.get(icre), "this contrib is zero since creation op destroys the bra");
        // if this is a particle moment contribution, the annihilation operator must be vacant in dst MBF
        DEBUG_ASSERT_TRUE(hole || !dst_mbf.get(icre), "this contrib is zero since annihilation op destroys the bra");

        auto& dst_finder = hole ? m_hole_dst_finder_table : m_particle_dst_finder_table;

        auto& dst_table = dst_finder.send(m_wf.m_dist.irank(dst_mbf));

        auto& send_row = dst_table.m_row;
        send_row.push_back_jump();

        send_row.m_dst_mbf = dst_mbf;
        send_row.m_hnci_factor = hnci_factor;
        send_row.m_ham_order = ham_order;
        send_row.m_ipart_dst = ipart_dst;
        send_row.m_ispinorb_left = ispinorb_left;
        send_row.m_ispinorb_right = ispinorb_right;
        return send_row;
    }

    /**
     * loop over all given spin-orbital indices and remove from mbf before calling the given functor
     * of course mbf will be modified each iteration, but we can promise constness since the mbf is as returned to
     * the initial state upon return
     */
    template<typename fn_t>
    static void foreach_occ(const Mbf& mbf, const uintv_t& ispinorbs, const fn_t& fn) {
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
    static void foreach_vac(const Mbf& mbf, const uintv_t& ispinorbs, const fn_t& fn) {
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

    /**
     * dispatch required version of the above loops
     */
    template<typename fn_t>
    static void foreach_perturber(const Mbf& mbf, const uintv_t& ispinorbs, bool hole, const fn_t& fn) {
        hole ? foreach_occ(mbf, ispinorbs, fn) : foreach_vac(mbf, ispinorbs, fn);
    }

public:
    SpecMoms(const Wavefunction& wf, Propagator& prop, const conf::SpecMoms& opts, sys::Basis basis):
        m_wf(wf), m_prop(prop), m_max_order(opts.m_max_order),
        m_hole_dst_finder_table("hole moment dst finder", DstFinderRow(basis), {1000, 1.0}),
        m_particle_dst_finder_table("particle moment dst finder", DstFinderRow(basis), {1000, 1.0}),
        m_selected_spinorbs(make_selected_spinorbs(opts, basis.m_frm.m_nspinorb)){
        for (uint_t ih=0ul; ih<m_max_order; ++ih) {
            m_hole_spec_moms.emplace_back(m_selected_spinorbs, basis.m_frm.m_nspinorb);
            m_particle_spec_moms.emplace_back(m_selected_spinorbs, basis.m_frm.m_nspinorb);
        }
    }

    void send_nondeterm_contribs(ExactLinear& prop, const field::FrmOnv& onv, bool hole, uint_t ispinorb_right) {
        conn::FrmOnv conn(onv);
        auto right_fn = [&](uint_t ispinorb_right) -> void {
//            auto mult_h_fn = [&]() -> void {
//
//            };
//            prop.conn_iters().loop(conn, onv, mult_h_fn);
        };
        foreach_perturber(onv, m_selected_spinorbs, hole, right_fn);
    }

    void send_nondeterm_contribs(
        StochLinear& prop,
        const field::FrmOnv& onv,
        wf_t weight,
        uint_t ipart_dst,
        bool hole,
        uint_t ispinorb_right)
    {
        uintv_t valid_spinorbs;
        valid_spinorbs.reserve(m_selected_spinorbs.size());
        for (auto ispinorb: m_selected_spinorbs) if (onv.get(ispinorb) != hole) valid_spinorbs.push_back(ispinorb);
        wf_t hnci_fac = 1.0;
        auto nattempt = prop.get_nattempt(weight);
        for (uint iattempt=0ul; iattempt < nattempt; ++iattempt) {

        }

        auto& exgens = prop.excit_gen_group();

        conn::FrmOnv conn(onv);
        auto right_fn = [&](uint_t ispinorb_right) -> void {
//            auto mult_h_fn = [&]() -> void {
//
//            };
//            prop.conn_iters().loop(conn, onv, mult_h_fn);
        };
        foreach_perturber(onv, m_selected_spinorbs, hole, right_fn);
    }

    void send_nondeterm_contribs(ExactLinear& prop, const field::FrmOnv& onv, bool hole) {
        conn::FrmOnv conn(onv);
        auto right_fn = [&](uint_t ispinorb_right) -> void {
//            auto mult_h_fn = [&]() -> void {
//
//            };
//            prop.conn_iters().loop(conn, onv, mult_h_fn);
        };
        foreach_perturber(onv, m_selected_spinorbs, hole, right_fn);
    }

    void send_nondeterm_contribs(StochLinear& prop, const field::FrmOnv& onv) {

    }

    void send_nondeterm_contribs(const field::FrmOnv& onv) {
        // dispatch correct overload based on the stochastic or exact type of the walker propagator
        {
            auto cast = dynamic_cast<ExactLinear*>(&m_prop);
            if (cast) send_nondeterm_contribs(*cast, onv);
        }
        {
            auto cast = dynamic_cast<StochLinear*>(&m_prop);
            if (cast) send_nondeterm_contribs(*cast, onv);
        }
    }
    void send_nondeterm_contribs(const field::BosOnv&) {}
    void send_nondeterm_contribs(const field::FrmBosOnv &onv) {send_nondeterm_contribs(onv.m_frm);}

    void make_nondeterm_contribs(const Table<DstFinderRow>& recv, v_t<SpecMom>& spec_moms) {
        auto dst_row = m_wf.m_store.m_row;
        for (auto row=recv.m_row; row; ++row) {
            m_wf.m_store.lookup(row.m_dst_mbf, dst_row);
            // check if the generated destination exists
            if (!dst_row) continue;
            const auto& dst_weight = dst_row.m_weight[row.m_ipart_dst];
            // or if the population on the replica is zero
            if (fptol::near_zero(dst_weight)) continue;
            // get the spectral moment order we are contributing to here
            auto& spec_mom = spec_moms[row.m_ham_order];
            spec_mom.make_contrib(row.m_ispinorb_left, row.m_ispinorb_right, dst_weight*row.m_hnci_factor);
        }
    }

    void make_nondeterm_contribs() {
        m_hole_dst_finder_table.communicate();
        make_nondeterm_contribs(m_hole_dst_finder_table.recv(), m_hole_spec_moms);
        m_particle_dst_finder_table.communicate();
        make_nondeterm_contribs(m_particle_dst_finder_table.recv(), m_particle_spec_moms);
    }


    operator bool() const {
        return !m_hole_spec_moms.empty() || !m_particle_spec_moms.empty();
    }
};
#endif //M7_SPECTRALMOMENT_H
