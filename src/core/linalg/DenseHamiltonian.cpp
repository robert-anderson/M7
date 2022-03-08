//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/Enumerators.h>
#include "DenseHamiltonian.h"
#include "src/core/foreach/Foreach.h"

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h, bool force_general) {
    if (force_general) return make_pair_iterator(h);
    typedef std::unique_ptr<PairBase> unique_t;
    if (!h.m_bd.m_nmode) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        auto fn = [this, &h](const field::FrmOnv &bra, size_t ibra, const field::FrmOnv &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };
        if (dynamic_cast<const SpinFrmHam *>(h.m_frm.get())) {
            return unique_t(new Pair<frm::Spins>({h.m_bd.m_nsite, h.m_frm->m_ms2_restrict}, fn));
        } else if (h.m_frm->m_kramers_attrs.conserving()) {
            return unique_t(new Pair<frm::Ms2Conserve>({h.m_bd.m_nsite, h.nelec(), h.m_frm->m_ms2_restrict}, fn));
        }
        return unique_t(new Pair<frm::General>({h.m_bd.m_nsite, h.nelec()}, fn));
    }
    else if (!h.m_bd.m_nsite) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
    }

    else {
        /*
         * hamiltonian is expressed in terms of fermion and boson operators, or it is assumed to be for testing purposes
         */
        bos::GeneralOpen bos_foreach(h.m_bd.m_nmode, h.m_nboson_max);
        auto fn = [this, &h](const field::FrmBosOnv &bra, size_t ibra, const field::FrmBosOnv &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };

        if (dynamic_cast<const SpinFrmHam *>(h.m_frm.get())) {
            typedef frm_bos::BosGeneralOpen<frm::Spins> foreach_t;
            return std::unique_ptr<PairBase>(
                    new Pair<foreach_t>({{h.m_bd.m_nsite, h.m_frm->m_ms2_restrict}, bos_foreach}, fn));
        } else if (h.m_frm->m_kramers_attrs.conserving()) {
            typedef frm_bos::BosGeneralOpen<frm::Ms2Conserve> foreach_t;
            return std::unique_ptr<PairBase>(
                    new Pair<foreach_t>({{h.m_bd.m_nsite, h.nelec(), h.m_frm->m_ms2_restrict}, bos_foreach}, fn));
        }
    }
    /*
     * the iterator has not already been returned, so return the force_general case appropriate for the basis dimensions
     */
    return make_pair_iterator(h);
}

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h) {
    if (!h.m_bd.m_nmode) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        auto fn = [this, &h](const field::FrmOnv &bra, size_t ibra, const field::FrmOnv &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };
        return std::unique_ptr<PairBase>(new Pair<frm::General>({h.m_bd.m_nsite, h.nelec()}, fn));
    } else if (!h.m_bd.m_nsite) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        auto fn = [this, &h](const field::BosOnv &bra, size_t ibra, const field::BosOnv &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };
        if (h.m_bos->m_nboson)
            return std::unique_ptr<PairBase>(new Pair<bos::GeneralClosed>({h.m_bd.m_nmode, h.m_bos->m_nboson}, fn));
        else
            return std::unique_ptr<PairBase>(new Pair<bos::GeneralOpen>({h.m_bd.m_nmode, h.m_nboson_max}, fn));
    } else {
        bos::GeneralOpen bos_foreach(h.m_bd.m_nmode, h.m_nboson_max);
        auto fn = [this, &h](const field::FrmBosOnv &bra, size_t ibra, const field::FrmBosOnv &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };
        typedef frm_bos::BosGeneralOpen<frm::General> foreach_t;
        return std::unique_ptr<PairBase>(
                new Pair<foreach_t>({{h.m_bd.m_nsite, h.nelec()}, bos_foreach}, fn));
    }
    ABORT("pair iterator not assigned");
    return nullptr;
}

size_t DenseHamiltonian::nrow(const Hamiltonian &h, bool force_general) {
    auto ptr = make_pair_iterator(h, force_general);
    return ptr->nrow();
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, bool force_general) :
    dense::SquareMatrix<defs::ham_t>(nrow(h, force_general)) {
    auto ptr = make_pair_iterator(h, force_general);
    ptr->loop();
}

