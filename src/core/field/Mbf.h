//
// Created by anderson on 11/29/21.
//

#ifndef M7_MBF_H
#define M7_MBF_H

#include <src/core/config/FciqmcConfig.h>
#include "Fields.h"

/**
 * namespace for utility methods that apply to all MBF types
 */

namespace mbf {

    void set_from_def(field::FrmOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);

    void set_from_def(field::FrmBosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);

    void set_from_def(field::BosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);

};


#endif //M7_MBF_H
