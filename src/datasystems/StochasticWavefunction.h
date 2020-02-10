//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_STOCHASTICWAVEFUNCTION_H
#define M7_STOCHASTICWAVEFUNCTION_H

#include "DataSystem.h"
#include "../data/PerforableMappedDataTable.h"

class StochasticWavefunction : public DataSystem {
    PerforableMappedDataTable store;
    DataTable send_buffer;
    DataTable recv_buffer;
public:
    //StochasticWavefunction();

};


#endif //M7_STOCHASTICWAVEFUNCTION_H
