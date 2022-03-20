//
// Created by rja on 17/08/2020.
//

#include "KramersSectorOccupation.h"

#if 0
KramersSectorOccupation::KramersSectorOccupation(size_t nelec) :
        m_nelec(nelec), m_sum(nelec+1){}

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

void KramersSectorOccupation::add(const DeterminantElement &det, const defs::wf_t &weight) {
    /*
     * n = 4
     *            s   (n+s)/2
     * 11110000   4      4
     * 11101000   2      3
     * 11001100   0      2
     * 10001110  -2      1
     * 00001111  -4      0
     *
     * n = 3
     *            s   (n+s)/2
     * 111000     3      3
     * 110100     1      2
     * 100110    -1      1
     * 000111    -3      0
     */
    m_sum[(m_nelec+det.spin())/2]+=std::abs(weight);
}

KramersSectorOccupation::KramersSectorOccupation(const DeterminantElement &det) :KramersSectorOccupation(det.nsetbit()){}

#endif