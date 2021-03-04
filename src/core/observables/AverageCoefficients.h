//
// Created by rja on 03/03/2021.
//

#ifndef M7_AVERAGECOEFFICIENTS_H
#define M7_AVERAGECOEFFICIENTS_H

#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "MevTable.h"

struct AverageCoefficients : BufferedTable<MevRow<defs::wf_t>, true>{

    AverageCoefficients(std::string name, defs::inds ninds, size_t nvalue, size_t nbucket=100):
        BufferedTable<MevRow<defs::wf_t>, true>(
                name,{{ninds, nvalue}, nbucket}){}

};


#endif //M7_AVERAGECOEFFICIENTS_H
