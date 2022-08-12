//
// Created by rja on 14/07/22.
//

#include <M7_lib/hamiltonian/frm/J1J2FrmHam.h>
#include "FciIters.h"


uint_t FciIters::niter_single() const {
    if (!m_single) return 0ul;
    return m_single->m_niter;
}

FciIters FciIters::make(const Hamiltonian& h, sys::Particles particles, bool force_general) {
    using namespace mbf_foreach;
    const sys::Sector sector(h.m_basis, particles);

    if (force_general) return make_general(h, particles);
    if (!sector.m_bos){
        /*
         * hamiltonian is boson operator-free, so can work in determinants: a.k.a. FrmOnvs
         */
        if (h.m_frm.is<SpinModelFrmHam>() || h.m_frm.is<J1J2FrmHam>()) return {frm::Spins(sector.m_frm)};
        else if (h.m_frm.m_kramers_attrs.conserving()) return {frm::Ms2Conserve(sector.m_frm)};
        return {frm::General(sector.m_frm)};
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
            if (h.m_frm.is<SpinModelFrmHam>() || h.m_frm.is<J1J2FrmHam>())
                return {frm_bos::OpenProduct<frm::Spins>(sector)};
            else if (sector.m_frm.m_elecs.m_ms2.conserve())
                return {frm_bos::OpenProduct<frm::Ms2Conserve>(sector)};
        }
    }
    /*
     * the iterator has not already been returned, so make with the force_general case appropriate for the basis dimensions
     */
    return make_general(h, particles);
}

FciIters FciIters::make(const Hamiltonian& h) {
    return make(h, h.default_particles(), false);
}

FciIters FciIters::make_general(const Hamiltonian& h, sys::Particles particles) {
    using namespace mbf_foreach;
    const sys::Sector sector(h.m_basis, particles);
    if (!sector.m_bos) {
        /*
         * hamiltonian is boson operator-free, can work in determinants: a.k.a. FrmOnvs
         */
        return {frm::General(sector.m_frm)};
    } else if (!sector.m_frm) {
        /*
         * hamiltonian is fermion operator-free, can work in permanents: a.k.a. BosOnvs
         */
        if (particles.m_bos.conserve()) return {bos::GeneralClosed(sector.m_bos)};
        else return {bos::GeneralOpen(sector.m_bos)};
    } else {
        if (sector.m_bos.m_bosons.conserve()) return {frm_bos::ClosedProduct<frm::General>(sector)};
        else return {frm_bos::OpenProduct<frm::General>(sector)};
    }
    return {};
}