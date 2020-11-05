//
// Created by rja on 06/07/2020.
//

#ifndef M7_HAMILTONIANCONNECTIONENUMERATOR_H
#define M7_HAMILTONIANCONNECTIONENUMERATOR_H


#include "src/core/basis/DeterminantConnection.h"
#include "src/core/hamiltonian/Hamiltonian.h"
#include "Enumerator.h"
#include "CombinationEnumerator.h"

#if 0

class HamiltonianSingleConnectionEnumerator : public Enumerator<MatrixElement<defs::ham_t>> {
    const Hamiltonian &m_h;
    const DeterminantElement &m_fonv;
    OccupiedOrbitals m_occs;
    VacantOrbitals m_vacs;

    size_t m_occind = ~0ul;
    size_t m_vacind = ~0ul;
    defs::ham_comp_t m_eps;

    MatrixElement<defs::ham_t> default_result() {
        return MatrixElement<defs::ham_t>(m_fonv);
    }

    bool next_element(MatrixElement<defs::ham_t> &result) override {
        if (++m_vacind<m_vacs.m_nind) {
            result.aconn.zero();
            result.aconn.add(
                    m_occs.m_inds[m_occind],
                    m_vacs.m_inds[m_vacind]
            );
            result.aconn.apply(m_fonv);
            result.element = m_h.get_element(result.aconn);
            if (consts::float_nearly_zero(result.element, m_eps)) {
                return next_element(result);
            }
        } else {
            m_vacind = ~0ul;
            if (++m_occind<m_occs.m_nind){
                return next_element(result);
            } else {
                m_occind = ~0ul;
                return false;
            }
        }
        return true;
    }

public:

    HamiltonianSingleConnectionEnumerator(const Hamiltonian &h, const DeterminantElement &det, const defs::ham_comp_t &eps=1e-12) :
            m_h(h), m_fonv(det), m_occs(det), m_vacs(det), m_eps(eps){
        m_occind++;
    }
};

class HamiltonianDoubleConnectionEnumerator : public Enumerator<MatrixElement<defs::ham_t>> {
    const Hamiltonian &m_h;
    const DeterminantElement &m_fonv;
    OccupiedOrbitals m_occs;
    VacantOrbitals m_vacs;

    CombinationEnumerator m_occ_enumerator;
    CombinationEnumerator m_vac_enumerator;

    defs::inds m_occinds, m_vacinds;
    defs::ham_comp_t m_eps;

    MatrixElement<defs::ham_t> default_result() {
        return MatrixElement<defs::ham_t>(m_fonv);
    }

    bool next_element(MatrixElement<defs::ham_t> &result) override {
        if (m_vac_enumerator.next(m_vacinds)) {
            result.aconn.zero();
            result.aconn.add(
                    m_occs.m_inds[m_occinds[0]],
                    m_occs.m_inds[m_occinds[1]],
                    m_vacs.m_inds[m_vacinds[0]],
                    m_vacs.m_inds[m_vacinds[1]]
            );
            result.element = m_h.get_element(result.aconn);
            if (consts::float_nearly_zero(result.element, m_eps)) {
                return next_element(result);
            }
        } else {
            if (m_occ_enumerator.next(m_occinds)) {
                return next_element(result);
            } else {
                return false;
            }
        }
        return true;
    }

public:

    HamiltonianDoubleConnectionEnumerator(const Hamiltonian &h, const DeterminantElement &det, const defs::ham_comp_t &eps=1e-12) :
            m_h(h), m_fonv(det), m_occs(det), m_vacs(det),
            m_occ_enumerator(m_occs.m_nind, 2),
            m_vac_enumerator(m_occs.m_nind, 2),
            m_occinds(2, ~0ul), m_vacinds(2, ~0ul), m_eps(eps){
        m_occ_enumerator.next(m_occinds);
    }
};


class HamiltonianConnectionEnumerator : public HamiltonianSingleConnectionEnumerator {
    HamiltonianDoubleConnectionEnumerator doubles_enumerator;
public:
    HamiltonianConnectionEnumerator(const Hamiltonian &h, const DeterminantElement &det, const defs::ham_comp_t &eps=1e-12) :
        HamiltonianSingleConnectionEnumerator(h, det), doubles_enumerator(h, det){
        m_subsequent = &doubles_enumerator;
    }
};

#if 0
AntisymConnection connection(ref);

FermionOnv excited(m_nsite);
for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
    const auto &occ = occs.m_inds[iocc];
    for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
        const auto &vac = vacs.m_inds[ivac];
        connection.zero();
        connection.add(occ, vac);
        connection.apply(ref, excited);
        auto helement = get_element(connection);
        if (!consts::float_nearly_zero(std::abs(helement), eps)) {
            size_t irow = list->push(excited);
            list->helement(irow) = helement;
        }
    }
}

ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, occs.m_nind, 2);
defs::inds occ_inds(2);
while (occ_enumerator.next(occ_inds)) {
    {
        ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
        defs::inds vac_inds(2);
        while (vac_enumerator.next(vac_inds)) {
            connection.zero();
            connection.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
            connection.apply(ref, excited);
            auto helement = get_element(connection);
            if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                size_t irow = list->push(excited);
                list->helement(irow) = helement;
                ASSERT(list->lookup(list->determinant(irow)) == irow);
            }
        }
    }
}
}
#endif

#endif //M7_HAMILTONIANCONNECTIONENUMERATOR_H
#endif //M7_HAMILTONIANCONNECTIONENUMERATOR_H
