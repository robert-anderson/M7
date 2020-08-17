//
// Created by rja on 17/08/2020.
//

#include "KramersSectorOccupation.h"

KramersSectorOccupation::KramersSectorOccupation(size_t nelec, size_t nsite) :
        m_nelec(nelec), m_nsite(nsite), m_sum(DeterminantElement::get_nspin_case(nelec)){}

KramersSectorOccupation::~KramersSectorOccupation() {
    std::cout << "Kramers sector occupation breakdown:" << std::endl;
    int mk = -m_nelec;
    for (auto& reducible: m_sum) {
        reducible.mpi_sum();
        if (mpi::i_am_root()) {
            std::cout << mk << " " << reducible.reduced() << std::endl;
            mk += 2;
        }
    }
}

size_t KramersSectorOccupation::add(const DeterminantElement &det, const defs::wf_t &weight) {
    m_sum[det.nsetbit()+det.spin()]+=std::abs(weight);
}

KramersSectorOccupation::KramersSectorOccupation(const DeterminantElement &det) :KramersSectorOccupation(det.nsetbit(), det.nsite()){}
