//
// Created by rja on 16/02/23.
//

#ifndef M7_SPECMOM_H
#define M7_SPECMOM_H


#include "M7_lib/util/Integer.h"
#include "M7_lib/field/Fields.h"
#include "M7_lib/propagator/Propagators.h"
#include "M7_lib/communication/Distribution.h"
#include "M7_lib/communication/SendRecv.h"

struct SpecMom {
    const uint_t m_order;
    /**
     * spin-orbital indices indexing the perturbers (a subset of all spinorbs)
     */
    const uintv_t m_selected_spinorbs;
    /**
     * the sampled quantities only have a memory requirement of creation values, so dense storage is fine
     */
    dense::Matrix<ham_t> m_store;
    /**
     * nspinorb element vector mapping perturber spinorbs gaplessly on [0, number of selected spinorbs)
     */
    const uintv_t m_spinorb_to_store_ind;

private:
    static uintv_t make_spinorb_to_store_ind(const uintv_t& selected_spinorbs, uint_t nspinorb) {
        uintv_t tmp(nspinorb, ~0ul);
        uint_t i = 0ul;
        for (auto ispinorb: selected_spinorbs) tmp[ispinorb] = i++;
        return tmp;
    }

public:

    SpecMom(uint_t order, uintv_t selected_spinorbs, uint_t nspinorb) :
        m_order(order), m_selected_spinorbs(std::move(selected_spinorbs)),
        m_store(m_selected_spinorbs.size(), m_selected_spinorbs.size()),
        m_spinorb_to_store_ind(make_spinorb_to_store_ind(m_selected_spinorbs, nspinorb)) {}

    void make_contrib(uint_t ispinorb_left, uint_t ispinorb_right, ham_t contrib) {
        const auto ileft = m_spinorb_to_store_ind[ispinorb_left];
        DEBUG_ASSERT_LT(ileft, m_store.nrow(), "left spin-orbital is not a selected perturber");
        const auto iright = m_spinorb_to_store_ind[ispinorb_right];
        DEBUG_ASSERT_LT(iright, m_store.nrow(), "right spin-orbital is not a selected perturber");
        m_store(ileft, iright) += contrib;
    }

    void save(hdf5::NodeWriter& gw) const {
        auto averaged = m_store;
        averaged.all_sum();
        averaged.save(convert::to_string(m_order), gw);
    }

};

/**
 * Row type is generic to hole/particle kinds. When a connection is generated, we must communicate the contribution
 * to the MPI rank that stores the instantaneous value of the destination MBF
 */
struct DstFinderRow : public Row {
    field::Mbf m_dst_mbf;
    field::Number<wf_t> m_hnci;
    field::Number<uint8_t> m_ham_order;
    field::Number<uint8_t> m_ipart_dst;
    field::Number<uint8_t> m_ispinorb_left;
    field::Number<uint8_t> m_ispinorb_right;

    explicit DstFinderRow(const sys::Basis& basis) :
            m_dst_mbf(this, basis, "destination MBF"),
            m_hnci(this, "source weight multiplied by H matrix elements and with the correct probabilistic normalization"),
            m_ham_order(this, "power of the Hamiltonian in the estimated expectation value"),
            m_ipart_dst(this, "WF part index of destination"),
            m_ispinorb_left(this, "index of the perturbing spin orbital at the left of the moment expectation value"),
            m_ispinorb_right(this, "index of the perturbing spin orbital at the right of the moment expectation value") {}
};

struct Comms {
    typedef SendRecv<DstFinderRow, Table<DstFinderRow>> dst_finder_comm_t;
    /**
     * communication tables for finding the dst weights in hole and particle-type moments
     */
    dst_finder_comm_t m_hole_dst_finder_table;
    dst_finder_comm_t m_particle_dst_finder_table;

    const Distribution& m_dist;

    Comms(sys::Basis basis, const Distribution& dist):
            m_hole_dst_finder_table("hole moment dst finder", DstFinderRow(basis), {1000, 1.0}),
            m_particle_dst_finder_table("particle moment dst finder", DstFinderRow(basis), {1000, 1.0}),
            m_dist(dist){}

    DstFinderRow& add_to_send(
            const field::Mbf& dst_mbf,
            ham_t hnci,
            uint_t ham_order,
            uint_t ipart_dst,
            uint_t ispinorb_left,
            uint_t ispinorb_right,
            bool hole)
    {
        auto& dst_finder = hole ? m_hole_dst_finder_table : m_particle_dst_finder_table;
        const auto dst_irank = m_dist.irank(dst_mbf);
        auto& dst_table = dst_finder.send(dst_irank);
        auto& send_row = dst_table.m_row;
        // if this is a hole moment contribution, the creation operator must be occupied in dst MBF
        DEBUG_ASSERT_FALSE(hole && !mbf::get_spinorb(dst_mbf, ispinorb_left),
                           "this contrib is zero since creation op destroys the bra");
        // if this is a particle moment contribution, the annihilation operator must be vacant in dst MBF
        DEBUG_ASSERT_FALSE(!hole && mbf::get_spinorb(dst_mbf, ispinorb_left),
                           "this contrib is zero since creation op destroys the bra");
        send_row.push_back_jump();

        send_row.m_dst_mbf = dst_mbf;
        send_row.m_hnci = hnci;
        send_row.m_ham_order = ham_order;
        send_row.m_ipart_dst = ipart_dst;
        send_row.m_ispinorb_left = ispinorb_left;
        send_row.m_ispinorb_right = ispinorb_right;
        return send_row;
    }

    void communicate() {
        m_hole_dst_finder_table.communicate();
        m_particle_dst_finder_table.communicate();
    }
};


struct Generator {
    const uintv_t& m_selected_spinorbs;
    Comms& m_comms;
    Generator(const uintv_t& selected_spinorbs, Comms& comms):
            m_selected_spinorbs(selected_spinorbs), m_comms(comms){}

    virtual void generate(const field::Mbf& src_mbf, wf_t weight, uint_t ipart_dst) = 0;
};

struct ExactGenerator : Generator {
    ExactLinear& m_prop;
    ExactGenerator(ExactLinear& prop, const uintv_t& selected_spinorbs, Comms& comms):
            Generator(selected_spinorbs, comms),
            m_prop(prop){}
};

struct StochGenerator : Generator {
    StochLinear& m_prop;
    PRNG m_prng;
    buffered::Mbf m_work_mbf;
    StochGenerator(StochLinear& prop, const uintv_t& selected_spinorbs, Comms& comms):
            Generator(selected_spinorbs, comms),
            m_prop(prop), m_prng(1, 100000), m_work_mbf(prop.m_ham.m_basis){}

    uint_t draw_ispinright() {
        return m_prng.draw_uint(m_selected_spinorbs.size());
    }

    void generate(const Mbf& src_mbf, wf_t weight, uint_t ipart_dst) override {
        const auto nattampt = m_prop.get_nattempt(weight);
        if (!nattampt) return;
        m_work_mbf = src_mbf;
        for (uint_t iattempt=0ul; iattempt < nattampt; ++iattempt) {
            auto ispinorb_right = draw_ispinright();
            (void) ispinorb_right;
            (void) ipart_dst;
            //auto ispinorb_leright = draw_ispinright();
        }
    }
};


#endif //M7_SPECMOM_H
