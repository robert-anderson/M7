//
// Created by Robert John Anderson on 2020-01-18.
//

#include <M7_lib/foreach/Foreach.h>

#include "DenseHamiltonian.h"

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h, bool force_general) {
    const auto& frm_hs = h.m_hs.m_frm;
    const auto& bos_hs = h.m_hs.m_bos;
    const auto nboson_conserve = bos_hs.m_nboson_conserve;

    if (force_general) return make_pair_iterator(h);
    if (!bos_hs) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        if (h.m_frm.is<SpinModelFrmHam>()){
            return unique_t(new frm::Pair<frm::Spins>(frm_hs));
        } else if (frm_hs.ms2_conserved()){
            return unique_t(new frm::Pair<frm::Ms2Conserve>(frm_hs));
        }
        return unique_t(new frm::Pair<frm::General>(frm_hs));
    } else if (!frm_hs) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
    } else {
        /*
         * hamiltonian is expressed in terms of fermion and boson operators, or it is assumed to be for testing purposes
         */
        if (nboson_conserve) {
            // closed system in the boson sector
        } else {
            // open system in the boson sector
            if (h.m_frm.is<SpinModelFrmHam>()) {
                return unique_t(new frm::Pair<frm::Spins>(frm_hs));
            }

            if (h.m_frm.is<SpinModelFrmHam>()) {
                typedef frm_bos::OpenProduct<frm::Spins> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                single_t foreach({frm_hs}, bos_hs);
                return unique_t(new pair_t(foreach));
            }
            else if (frm_hs.ms2_conserved()) {
                typedef frm_bos::OpenProduct<frm::Ms2Conserve> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                single_t foreach({frm_hs}, bos_hs);
                return unique_t(new pair_t(foreach));
            }
        }
    }
    /*
     * the iterator has not already been returned, so make with the force_general case appropriate for the basis dimensions
     */
    return make_pair_iterator(h);
}

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h) {
    const auto& frm_hs = h.m_hs.m_frm;
    const auto& bos_hs = h.m_hs.m_bos;
    const auto nboson_conserve = bos_hs.m_nboson_conserve;
    if (!bos_hs) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        return unique_t(new frm::Pair<frm::General>(frm_hs));
    } else if (!frm_hs) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        if (nboson_conserve)
            return unique_t(new bos::Pair<bos::GeneralClosed>(bos_hs));
        else
            return unique_t(new bos::Pair<bos::GeneralOpen>(bos_hs));
    } else {
        if (nboson_conserve) {
            typedef frm_bos::ClosedProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            single_t foreach({frm_hs}, bos_hs);
            return unique_t(new pair_t(foreach));
        }
        else {
            typedef frm_bos::OpenProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            single_t foreach({frm_hs}, bos_hs);
            return unique_t(new pair_t(foreach));
        }
    }
    ABORT("pair iterator not assigned");
}

size_t DenseHamiltonian::nrow(const Hamiltonian &h, bool force_general) {
    auto ptr = make_pair_iterator(h, force_general);
    return ptr->m_nrow;
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, bool force_general) :
        dense::SquareMatrix<defs::ham_t>(nrow(h, force_general)) {
    auto ptr = make_pair_iterator(h, force_general);
    /*
     * only one of the virtual methods will be overridden to a non-empty loop so they can all be called
     */
    loop_over_pair_iterator<field::FrmOnv>(ptr.get(), h);
    loop_over_pair_iterator<field::BosOnv>(ptr.get(), h);
    loop_over_pair_iterator<field::FrmBosOnv>(ptr.get(), h);
}