//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>

class FciqmcCalculation {

public:
    FciqmcCalculation(){

        auto r = omp_get_num_threads();


    }

    void execute();

};


#endif //M7_FCIQMCCALCULATION_H
