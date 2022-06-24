//
// Created by Robert John Anderson on 2020-01-18.
//

#include "DenseHamiltonian.h"

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(
        const Hamiltonian &h, sys::Particles particles, bool force_general) {
    const sys::Sector sector(h.m_basis, particles);

    if (force_general) return make_pair_iterator(h, particles);
    if (!sector.m_bos){
        /*
         * hamiltonian is boson operator-free, so can work in determinants: a.k.a. FrmOnvs
         */
        if (h.m_frm.is<SpinModelFrmHam>()){
            return unique_t(new frm::Pair<frm::Spins>(sector.m_frm));
        } else if (h.m_frm.m_kramers_attrs.conserving()){
            return unique_t(new frm::Pair<frm::Ms2Conserve>(sector.m_frm));
        }
        return unique_t(new frm::Pair<frm::General>(sector.m_frm));
    } else if (!sector.m_frm) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
    } else {
        /*
         * hamiltonian is expressed in terms of fermion and boson operators, or it is assumed to be for testing purposes
         */
        if (h.m_boson_number_conserve) {
            // closed system in the boson sector
        } else {
            // open system in the boson sector
            if (h.m_frm.is<SpinModelFrmHam>()) {
                typedef frm_bos::OpenProduct<frm::Spins> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                return unique_t(new pair_t(sector));
            }
            else if (sector.m_frm.m_elecs.m_ms2.conserve()) {
                typedef frm_bos::OpenProduct<frm::Ms2Conserve> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                return unique_t(new pair_t(sector));
            }
        }
    }
    /*
     * the iterator has not already been returned, so make with the force_general case appropriate for the basis dimensions
     */
    return make_pair_iterator(h, particles);
}

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h, sys::Particles particles) {

    const sys::Sector sector(h.m_basis, particles);
    if (!sector.m_bos) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        return unique_t(new frm::Pair<frm::General>(sector.m_frm));
    } else if (!sector.m_frm) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        if (particles.m_bos.conserve())
            return unique_t(new bos::Pair<bos::GeneralClosed>(sector.m_bos));
        else
            return unique_t(new bos::Pair<bos::GeneralOpen>(sector.m_bos));
    } else {
        if (sector.m_bos.m_bosons.conserve()) {
            typedef frm_bos::ClosedProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            return unique_t(new pair_t(sector));
        }
        else {
            typedef frm_bos::OpenProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            return unique_t(new pair_t(sector));
        }
    }
    ABORT("pair iterator not assigned");
    return nullptr;
}

uint_t DenseHamiltonian::nrow(const Hamiltonian &h, sys::Particles particles, bool force_general) {
    auto ptr = make_pair_iterator(h, particles, force_general);
    return ptr->m_nrow;
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, sys::Particles particles, bool force_general) :
        dense::SquareMatrix<ham_t>(nrow(h, particles, force_general)) {
    auto ptr = make_pair_iterator(h, particles, force_general);

    if (h.m_basis.m_frm && h.m_basis.m_bos) {
        buffered::FrmBosOnv bra(h.m_basis);
        auto ket = bra;
        loop_over_pair_iterator(ptr.get(), h, bra, ket);
    }
    else if (h.m_basis.m_frm) {
        buffered::FrmOnv bra(h.m_basis.m_frm);
        auto ket = bra;
        loop_over_pair_iterator(ptr.get(), h, bra, ket);
    }
    else if (h.m_basis.m_bos) {
        buffered::BosOnv bra(h.m_basis.m_bos);
        auto ket = bra;
        loop_over_pair_iterator(ptr.get(), h, bra, ket);
    }
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, bool force_general) :
        DenseHamiltonian(h, h.default_particles(), force_general){}
