//
// Created by anderson on 11/29/21.
//

#ifndef M7_MBF_H
#define M7_MBF_H

#include <M7_lib/config/FciqmcConfig.h>

#include "Fields.h"

/**
 * namespace for utility methods that apply to all MBF types
 */

namespace mbf {

    /**
     * set the referenced ONV object to the assumed Hartree--Fock determinant within the given sector
     */
    void set_aufbau_mbf(field::FrmOnv &onv);

    /**
     * there is no proper aufbau principle for bosons, so fill from the 0 mode until max nbosons is reached, then
     * advance to the next mode
     */
    void set_aufbau_mbf(field::BosOnv &onv);


    /**
     * set the referenced ONV object to the "anti-ferromagnetic" configuration
     */
    void set_neel_mbf(field::FrmOnv &onv);

    void set_from_def_array(field::FrmOnv &mbf, const std::vector<defs::inds> &def, size_t idef);

    void set_from_def_array(field::BosOnv &mbf, const std::vector<defs::inds> &def, size_t idef);

    void set(field::FrmOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);

    void set(field::BosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);

    void set(field::FrmBosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef);
};


#endif //M7_MBF_H
