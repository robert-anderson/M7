//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"

void FciqmcCalculation::run() {
    //for (auto i{0ul}; i<1000; ++i) {
    while(1){
        m_wf.evolve(m_p);
    }
}
