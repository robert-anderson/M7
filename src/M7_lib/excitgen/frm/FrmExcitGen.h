//
// Created by Robert J. Anderson on 03/04/2022.
//

#ifndef M7_FRMEXCITGEN_H
#define M7_FRMEXCITGEN_H

#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct FrmExcitGen : ExcitGen {

    const FrmHam& m_h;

    FrmExcitGen(const FrmHam& h, PRNG& prng, uintv_t exsigs, str_t description);

    bool draw_frmbos(uint_t exsig, const field::FrmBosOnv& src,
                     prob_t& prob, conn::FrmBosOnv& conn) override;

    bool draw_h_frm(uint_t exsig, const field::FrmOnv& src, prob_t& prob,
                    ham_t& helem, conn::FrmOnv& conn) override;

    bool draw_h_frmbos(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob,
                       ham_t& helem, conn::FrmBosOnv& conn) override;

    bool draw_h_bos(uint_t /*exsig*/, const field::BosOnv& /*src*/, prob_t& prob,
                    ham_t& helem, conn::BosOnv& /*conn*/) override {
        prob = 0.0;
        helem = 0.0;
        return false;
    }

};

struct FrmLatticeExcitGen : FrmExcitGen {
protected:
    mutable lattice::adj_row_t m_work_adj_row;
    mutable v_t<const lattice::AdjElement*> m_valid_in_adj_row;

    /**
     * set the vector of pointers to "valid" elements of the adjacent row. the validity criterion is determined by the
     * logic prescribed in the is_valid_fn functor
     * @tparam is_valid_fn_t
     *  functor which returns true if the given site is to be included in the list of valid elements of m_work_adj_row
     * @param is_valid_fn
     *  instance of the functor
     */
    template<typename is_valid_fn_t>
    void set_valid_adj(const is_valid_fn_t& is_valid_fn) const {
        functor::assert_prototype<bool(uint_t)>(is_valid_fn);
        m_valid_in_adj_row.clear();
        for (const auto& adj : m_work_adj_row) {
            if (is_valid_fn(adj.m_isite)) m_valid_in_adj_row.push_back(&adj);
        }
    }
    /**
     * convenient wrapper for the above general method, which includes those elements of m_work_adj_row which represent
     * indices of sites which are unoccupied the src determinant in the given spin channel
     * @param src
     *  fermion ONV
     * @param ispin
     *  spin channel in which to look for occupied spin orbitals
     */
    void set_valid_adj_vacant(const field::FrmOnv& src, uint_t ispin) const;


public:
    FrmLatticeExcitGen(const FrmHam& h, PRNG& prng, uintv_t exsigs, str_t description):
        FrmExcitGen(h, prng, exsigs, description){
        REQUIRE_TRUE(m_h.m_basis.m_lattice.get(), "Lattice excitation generator requires lattice definition in basis");
        const auto nadj_max = m_h.m_basis.m_lattice->m_nadj_max;
        m_work_adj_row.reserve(nadj_max);
        m_valid_in_adj_row.reserve(nadj_max);
    }
};


#endif //M7_FRMEXCITGEN_H
