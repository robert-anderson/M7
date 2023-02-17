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

    void set_from_def_array(field::FrmOnv &mbf, const v_t<uintv_t> &def, uint_t idef);

    void set_from_def_array(field::BosOnv &mbf, const v_t<uintv_t> &def, uint_t idef);

    void set(field::FrmOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef);

    void set(field::BosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef);

    void set(field::FrmBosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef);

    OpSig exsig(const field::FrmOnv &src, const field::FrmOnv &dst);

    OpSig exsig(const field::BosOnv &src, const field::BosOnv &dst);

    OpSig exsig(const field::FrmBosOnv &src, const field::FrmBosOnv &dst);

    static bool get_spinorb(const field::FrmOnv& onv, uint_t ispinorb) {
        return onv.get(ispinorb);
    }
    static bool get_spinorb(const field::BosOnv&, uint_t) {
        return false;
    }
    static bool get_spinorb(const field::FrmBosOnv& onv, uint_t ispinorb) {
        return get_spinorb(onv.m_frm, ispinorb);
    }

    static void put_spinorb(field::FrmOnv& onv, uint_t ispinorb, bool v) {
        onv.put(ispinorb, v);
    }
    static void put_spinorb(field::BosOnv&, uint_t, bool){}
    static void put_spinorb(field::FrmBosOnv& onv, uint_t ispinorb, bool v) {
        put_spinorb(onv.m_frm, ispinorb, v);
    }
};


#endif //M7_MBF_H
