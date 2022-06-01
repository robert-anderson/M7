//
// Created by Robert John Anderson on 2020-01-18.
//

#include <M7_lib/foreach/Foreach.h>

#include "DenseHamiltonian.h"

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(
        const Hamiltonian &h, sys::Particles particles, bool force_general) {
    const sys::frm::Sector frm_sector(h.m_basis.m_frm, particles.m_frm);
    const sys::bos::Sector bos_sector(h.m_basis.m_bos, particles.m_bos);

    if (force_general) return make_pair_iterator(h, particles);
    if (!bos_sector){
        /*
         * hamiltonian is boson operator-free, so can work in determinants: a.k.a. FrmOnvs
         */
        if (h.m_frm.is<SpinModelFrmHam>()){
            return unique_t(new frm::Pair<frm::Spins>(frm_sector));
        } else if (h.m_frm.m_kramers_attrs.conserving()){
            return unique_t(new frm::Pair<frm::Ms2Conserve>(frm_sector));
        }
        return unique_t(new frm::Pair<frm::General>(frm_sector));
    } else if (!frm_sector) {
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
                return unique_t(new frm::Pair<frm::Spins>(frm_sector));
            }

            if (h.m_frm.is<SpinModelFrmHam>()) {
                typedef frm_bos::OpenProduct<frm::Spins> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                single_t foreach({frm_sector}, bos_sector);
                return unique_t(new pair_t(foreach));
            }
            else if (frm_sector.m_elecs.m_ms2.conserve()) {
                typedef frm_bos::OpenProduct<frm::Ms2Conserve> single_t;
                typedef frm_bos::Pair<single_t> pair_t;
                single_t foreach({frm_sector}, bos_sector);
                return unique_t(new pair_t(foreach));
            }
        }
    }
    /*
     * the iterator has not already been returned, so make with the force_general case appropriate for the basis dimensions
     */
    return make_pair_iterator(h, particles);
}

std::unique_ptr<PairBase> DenseHamiltonian::make_pair_iterator(const Hamiltonian &h, sys::Particles particles) {
    const sys::frm::Sector frm_sector(h.m_basis.m_frm, particles.m_frm);
    const sys::bos::Sector bos_sector(h.m_basis.m_bos, particles.m_bos);

    if (!h.m_bos) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        return unique_t(new frm::Pair<frm::General>(frm_sector));
    } else if (!h.m_frm) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        if (particles.m_bos.conserve())
            return unique_t(new bos::Pair<bos::GeneralClosed>(bos_sector));
        else
            return unique_t(new bos::Pair<bos::GeneralOpen>(bos_sector));
    } else {
        if (bos_sector.m_bosons.conserve()) {
            typedef frm_bos::ClosedProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            single_t foreach({frm_sector}, bos_sector);
            return unique_t(new pair_t(foreach));
        }
        else {
            typedef frm_bos::OpenProduct<frm::General> single_t;
            typedef frm_bos::Pair<single_t> pair_t;
            single_t foreach({frm_sector}, bos_sector);
            return unique_t(new pair_t(foreach));
        }
    }
    ABORT("pair iterator not assigned");
    return nullptr;
}

size_t DenseHamiltonian::nrow(const Hamiltonian &h, sys::Particles particles, bool force_general) {
    auto ptr = make_pair_iterator(h, particles, force_general);
    return ptr->m_nrow;
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, sys::Particles particles, bool force_general) :
        dense::SquareMatrix<defs::ham_t>(nrow(h, particles, force_general)) {
    auto ptr = make_pair_iterator(h, particles, force_general);
    /*
     * only one of the virtual methods will be overridden to a non-empty loop so they can all be called
     */
    loop_over_pair_iterator<field::FrmOnv>(ptr.get(), h);
    loop_over_pair_iterator<field::BosOnv>(ptr.get(), h);
    loop_over_pair_iterator<field::FrmBosOnv>(ptr.get(), h);
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, bool force_general) :
        DenseHamiltonian(h, h.default_particles(), force_general){}
