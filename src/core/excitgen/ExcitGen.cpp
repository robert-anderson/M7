//
// Created by rja on 04/06/2020.
//

#include "ExcitGen.h"

FrmExcitGen::FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t exsig) :
ExcitGen(h, prng, {exsig}),
m_spin_conserving(exsig == exsig_utils::ex_single ? h.m_frm.m_kramers_attrs.m_conserving_singles : h.m_frm.m_kramers_attrs.m_conserving_double) {}

bool FrmExcitGen::draw(const FrmOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                       defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
    return ExcitGen::draw(src_onv, occs, vacs, prob, helem, conn);
}

bool FrmExcitGen::draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                       defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
    return draw(src_onv.m_frm, occs, vacs, prob, helem, conn.m_frm);
}

ExcitGen::ExcitGen(const Hamiltonian &h, PRNG &prng, defs::inds exsigs) :
        m_h(h), m_prng(prng),
        m_nspinorb(m_h.nsite() * 2),
        m_nelec(m_h.nelec()),
        m_norb_pair(integer_utils::nspair(m_nspinorb)),
        m_nelec_pair(integer_utils::nspair(m_nelec)), m_exsigs(std::move(exsigs)) {
}

bool ExcitGen::draw(const FrmOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
                    defs::ham_t &helem, conn::FrmOnv &conn) {
    prob = 0.0;
    helem = 0.0;
    return false;
}

bool ExcitGen::draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs, defs::prob_t &prob,
               defs::ham_t &helem, conn::FrmBosOnv &conn) {
    prob = 0.0;
    helem = 0.0;
    return false;
}
