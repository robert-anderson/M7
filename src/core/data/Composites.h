//
// Created by RJA on 01/11/2020.
//

#ifndef M7_COMPOSITES_H
#define M7_COMPOSITES_H

#include "ConfigurationField.h"

namespace composites {

    using Configuration = ConfigurationField<0ul>;

    template<size_t nind>
    using Configurations = ConfigurationField<nind>;
};


#endif //M7_COMPOSITES_H
