//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include "DataSystem.h"
#include "src/data/PerforableMappedTable.h"

class Wavefunction {
    PerforableMappedTable<Determinant> store;
    Table send_buffer;
    Table recv_buffer;
public:

    /*
    Wavefunction(const Determinant& reference){

    }*/
/*
    void propagate() {
#pragma omp parallel {
        Table thread_send_buffer;
#pragma omp for
        for (auto irow{0ul}; irow < store.highwatermark()[0]; ++irow) {

        }
    }
}*/

};

#endif //M7_WAVEFUNCTION_H
