//
// Created by Robert J. Anderson on 11/29/21.
//

#ifndef M7_MBF_H
#define M7_MBF_H

#include <M7_lib/conf/Conf.h>

#include "Fields.h"

/**
 * namespace for utility methods that apply to all MBF types
 */

namespace mbf {

    /**
     * set the referenced ONV object to the assumed Hartree--Fock determinant within the sector specified by the
     * elecs argument
     */
    void set_aufbau_mbf(field::FrmOnv &onv, sys::frm::Electrons elecs);
    void set_aufbau_mbf(field::FrmOnv &onv, sys::Particles particles);

    /**
     * there is no proper aufbau principle for bosons, so fill from the 0 mode until max nbosons is reached, then
     * advance to the next mode until the number of quanta specified by the bosons argument has been reached
     */
    void set_aufbau_mbf(field::BosOnv &onv, sys::bos::Bosons bosons);
    void set_aufbau_mbf(field::BosOnv &onv, sys::Particles particles);

    /**
     * set the referenced ONV object to the "anti-ferromagnetic" configuration
     */
    void set_neel_mbf(field::FrmOnv &onv, sys::frm::Electrons elecs);

    void set_from_def_array(field::FrmOnv &mbf, const std::vector<defs::inds> &def, size_t idef);

    void set_from_def_array(field::BosOnv &mbf, const std::vector<defs::inds> &def, size_t idef);

    void set(field::FrmOnv &mbf, sys::Particles particles, const conf::MbfDef &def, size_t idef);

    void set(field::BosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, size_t idef);

    void set(field::FrmBosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, size_t idef);
};


#endif //M7_MBF_H
