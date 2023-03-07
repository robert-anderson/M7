//
// Created by rja on 16/02/23.
//

#ifndef M7_SPECMOMS_H
#define M7_SPECMOMS_H

#include "SpecMom.h"
#include "M7_lib/util/Vector.h"

struct SpecMoms {

    const conf::SpecMoms& m_opts;

    typedef SendRecv<DstFinderRow, Table<DstFinderRow>> dst_finder_comm_t;
    /**
     * const pointer to wavefunction class so we can lookup dst walker weights
     */
    const wf::Vectors* m_wf = nullptr;
    /**
     * pointer to propagator class which is used to generate connections in the perturbed many-body basis
     */
    Propagator* m_prop = nullptr;
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

    const Epoch& m_accum_epoch;


//    DstFinderRow& add_send_contrib(
//            bool hole,
//            const Mbf& dst_mbf,
//            uint_t irank_dst,
//            ham_t hnci,
//            uint_t ham_order,
//            uint_t ipart_dst,
//            uint_t ispinorb_left,
//            uint_t ispinorb_right) {
//        // if this is a hole moment contribution, the creation operator must be occupied in dst MBF
//        DEBUG_ASSERT_FALSE(hole && !mbf::get_spinorb(dst_mbf, ispinorb_left),
//                           "this contrib is zero since creation op destroys the bra");
//        // if this is a particle moment contribution, the annihilation operator must be vacant in dst MBF
//        DEBUG_ASSERT_FALSE(!hole && mbf::get_spinorb(dst_mbf, ispinorb_left),
//                           "this contrib is zero since creation op destroys the bra");
//
//        auto& dst_finder = hole ? m_hole_dst_finder_table : m_particle_dst_finder_table;
//        auto& dst_table = dst_finder.send(irank_dst);
//
//        auto& send_row = dst_table.m_row;
//        send_row.push_back_jump();
//
//        send_row.m_dst_mbf = dst_mbf;
//        send_row.m_hnci = hnci;
//        send_row.m_ham_order = ham_order;
//        send_row.m_ipart_dst = ipart_dst;
//        send_row.m_ispinorb_left = ispinorb_left;
//        send_row.m_ispinorb_right = ispinorb_right;
//        return send_row;
//    }

private:
    static uintv_t make_selected_spinorbs(const conf::SpecMoms& opts, uint_t nspinorb) {
        if (opts.m_spinorbs.m_value.empty()) return vector::range(0ul, nspinorb);
        return opts.m_spinorbs;
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
    SpecMoms(const conf::SpecMoms& opts, sys::Sector sector, const Epoch& accum_epoch):
        m_opts(opts), m_max_order(opts.m_max_order),
        m_hole_dst_finder_table("hole moment dst finder", DstFinderRow(sector.basis()), {1000, 1.0}),
        m_particle_dst_finder_table("particle moment dst finder", DstFinderRow(sector.basis()), {1000, 1.0}),
        m_selected_spinorbs(make_selected_spinorbs(opts, sector.basis().m_frm.m_nspinorb)),
        m_accum_epoch(accum_epoch) {
        if (!m_opts.m_enabled) return;
        for (uint_t order = 0ul; order <= m_max_order; ++order) {
            m_hole_spec_moms.emplace_back(order, m_selected_spinorbs, sector.m_frm.m_basis.m_nspinorb);
            m_particle_spec_moms.emplace_back(order, m_selected_spinorbs, sector.m_frm.m_basis.m_nspinorb);
        }
    }

    ~SpecMoms() {
        if (m_opts.m_save.m_enabled) save();
    }

    operator bool() const {
        return !m_hole_spec_moms.empty() || !m_particle_spec_moms.empty();
    }

    void save(const hdf5::NodeWriter& parent) {
        if (!m_accum_epoch) {
            logging::warn("MAE accumulation epoch was not reached in this calculation: omitting spectral moment save");
            return;
        }
        hdf5::GroupWriter archive(parent, "archive");
        hdf5::GroupWriter hole_group(archive, "hole");
        hdf5::GroupWriter particle_group(archive, "particle");
        {
            // unnormalized SpecMoms, suitable for restarts
            for (const auto& specmom: m_hole_spec_moms) specmom.save(hole_group);
            for (const auto& specmom: m_particle_spec_moms) specmom.save(particle_group);
        }
    }

    void save() {
        if (!m_accum_epoch) {
            logging::warn("MAE accumulation epoch was not reached in this calculation: omitting SpecMoms save");
            return;
        }
        REQUIRE_TRUE(m_opts.m_save.m_enabled, "save() called on SpecMoms object but saving was not enabled");
        save(hdf5::FileWriter("M7.spec_mom.h5"));
    }

};


#endif //M7_SPECMOMS_H
