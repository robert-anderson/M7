//
// Created by Robert John Anderson on 2020-01-18.
//

#include <M7_lib/foreach/Foreach.h>

#include "DenseHamiltonian.h"

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h, bool force_general) {
    const auto nsite = h.m_bd.m_frm.m_nsite;
    const auto nelec = h.m_frm.m_hd.m_nelec;
    const auto ms2_restrict = h.m_frm.m_hd.m_ms2_restrict;
    const auto nmode = h.m_bd.m_bos.m_nmode;
    const auto nboson = h.m_bos.m_hd.m_nboson;
    const auto nboson_max = h.m_bos.m_hd.m_nboson_max;

    if (force_general) return make_pair_iterator(h);
    if (!nmode) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        if (dynamic_cast<const SpinModelFrmHam *>(&h.m_frm)) {
            return unique_t(new frm::Pair<frm::Spins>({nsite, ms2_restrict}));
        } else if (h.m_frm.m_kramers_attrs.conserving()) {
            return unique_t(new frm::Pair<frm::Ms2Conserve>({nsite, h.m_frm.m_hd}));
        }
        return unique_t(new frm::Pair<frm::General>({nsite, nelec}));
    } else if (!nsite) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
    } else {
        /*
         * hamiltonian is expressed in terms of fermion and boson operators, or it is assumed to be for testing purposes
         */
        if (!nboson) {
            // open system in the boson sector
            if (dynamic_cast<const SpinModelFrmHam *>(&h.m_frm)) {
                return unique_t(new frm::Pair<frm::Spins>({nsite, ms2_restrict}));
            }
        } else {
            // closed system in the boson sector

        }

        if (dynamic_cast<const SpinModelFrmHam *>(&h.m_frm)) {
            typedef frm_bos::OpenProduct<frm::Spins> foreach_t;
            return unique_t(new frm_bos::Pair<foreach_t>({{nsite, ms2_restrict}, nmode, nboson_max}));
        } else if (h.m_frm.m_kramers_attrs.conserving()) {
            typedef frm_bos::OpenProduct<frm::Ms2Conserve> foreach_t;
            return unique_t(new frm_bos::Pair<foreach_t>({{nsite, h.m_frm.m_hd}, nmode, nboson_max}));
        }
    }
    /*
     * the iterator has not already been returned, so make with the force_general case appropriate for the basis dimensions
     */
    return make_pair_iterator(h);
}

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h) {
    const auto nsite = h.m_bd.m_frm.m_nsite;
    const auto nelec = h.m_frm.m_hd.m_nelec;
    const auto nmode = h.m_bd.m_bos.m_nmode;
    const auto nboson = h.m_bos.m_hd.m_nboson;
    const auto nboson_max = h.m_bos.m_hd.m_nboson_max;
    if (!nmode) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        return unique_t(new frm::Pair<frm::General>({nsite, nelec}));
    } else if (!nsite) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        if (nboson)
            return unique_t(new bos::Pair<bos::GeneralClosed>({nmode, nboson}));
        else
            return unique_t(new bos::Pair<bos::GeneralOpen>({nmode, nboson_max}));
    } else {
        if (nboson) {
            typedef frm_bos::OpenProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            return unique_t(new pair_t({{nsite, nelec}, nmode, nboson}));
        } else {
            typedef frm_bos::ClosedProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            return unique_t(new pair_t({{nsite, nelec}, nmode, nboson_max}));
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

