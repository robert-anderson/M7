//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <M7_lib/conf/Conf.h>
#include <M7_lib/wavefunction/Wavefunction.h>

#include "M7_lib/propagator/StochLinear.h"
#include "Solver.h"

namespace fciqmc {
     /**
      * @param opts
      *  configuration document specifying the calculation to be performed
      */
    void run(const conf::Document& opts);
}

#endif //M7_FCIQMCCALCULATION_H