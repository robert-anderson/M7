//
// Created by rja on 05/07/2021.
//

#ifndef M7_ARCHIVABLE_H
#define M7_ARCHIVABLE_H

#include "HDF5Wrapper.h"


struct Archivable {
    const std::string m_name;
    /**
     * include this object in checkpoint saves
     */
    const bool m_load;
    /**
     * include this object in finalizing saves
     */
    const bool m_save;
    /**
     * include this object in checkpoint saves
     */
    const bool m_chkpt;

    Archivable(std::string name, bool load, bool save, bool chkpt):
        m_name(name), m_load(load), m_save(save), m_chkpt(chkpt){}

    virtual ~Archivable(){}

    virtual void load_fn(hdf5::GroupReader& parent) = 0;
    virtual void save_fn(hdf5::GroupWriter& parent) = 0;

    void load(hdf5::GroupReader& parent) {
        if(m_load) {
            REQUIRE_TRUE_ALL(parent.child_exists(m_name), log::format("Load failed: \"{}\" does not exist in archive", m_name));
            load_fn(parent);
        }
    }
    void save(hdf5::GroupWriter& parent){
        if (m_save) save_fn(parent);
    }
    void chkpt(hdf5::GroupWriter& parent){
        if (m_save) save_fn(parent);
    }

};


#endif //M7_ARCHIVABLE_H
