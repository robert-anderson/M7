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
private:
    typedef v_t<const lattice::Lattice::adj_t*> valid_adj_t;
    mutable valid_adj_t m_work_valid_adj;
protected:
    /**
     * set the vector of pointers to "valid" elements of the adjacent row. the validity criterion is determined by the
     * logic prescribed in the is_valid_fn functor
     * @tparam is_valid_fn_t
     *  functor which returns true if the given site is to be included in the list of valid elements of m_work_adj_row
     * @param is_valid_fn
     *  instance of the functor
     * @return
     *  const ref to filled vector of pointers to valid adjacent site indices
     */
    template<typename is_valid_fn_t>
    const valid_adj_t& valid_adj(uint_t isite, const is_valid_fn_t& is_valid_fn) const {
        functor::assert_prototype<bool(uint_t)>(is_valid_fn);
        m_work_valid_adj.clear();
        const auto& adj = m_h.m_basis.m_lattice->m_sparse_adj;
        for (auto it=adj.cbegin(isite); it!=adj.cend(isite); ++it){
            if (is_valid_fn(it->first)) m_work_valid_adj.push_back(it.base());
        }
        return m_work_valid_adj;
    }
    /**
     * convenient wrapper for the above general method, which includes those elements of m_work_adj_row which represent
     * indices of sites which are unoccupied the src determinant in the given spin channel
     * @param src
     *  fermion ONV
     * @param ispin
     *  spin channel in which to look for occupied spin orbitals
     * @return
     *  const ref to filled vector of pointers to valid adjacent site indices
     */
    const valid_adj_t& valid_adj(uint_t isite, const field::FrmOnv& src, uint_t ispin) const;


public:
    FrmLatticeExcitGen(const FrmHam& h, PRNG& prng, uintv_t exsigs, str_t description):
        FrmExcitGen(h, prng, exsigs, description){
        REQUIRE_TRUE(m_h.m_basis.m_lattice.get(), "Lattice excitation generator requires lattice definition in basis");
        const auto nadj_max = m_h.m_basis.m_lattice->m_nadj_max;
        m_work_valid_adj.reserve(nadj_max);
    }
};


#endif //M7_FRMEXCITGEN_H
