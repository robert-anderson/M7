//
// Created by rja on 05/07/2021.
//

#ifndef M7_ARCHIVABLE_H
#define M7_ARCHIVABLE_H

#include <M7_lib/config/FciqmcConfig.h>
#include <M7_lib/util/Timer.h>

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

    Archivable(std::string name, bool load, bool save, bool chkpt);

    Archivable(std::string name, const fciqmc_config::Archivable &opts);

    /**
     * ctor for object which are not optionally archivable
     * @param name
     *  label of group in HDF5 file
     * @param archive_opts
     *  Global options for the archive
     */
    Archivable(std::string name, const fciqmc_config::Archive &archive_opts);

    virtual ~Archivable() {}

protected:
    virtual void load_fn(hdf5::GroupReader &parent) = 0;

    virtual void save_fn(hdf5::GroupWriter &parent) = 0;

public:
    void load(hdf5::GroupReader &parent);

    void save(hdf5::GroupWriter &parent);

    void chkpt(hdf5::GroupWriter &parent);
};

struct Metadata : Archivable {
    size_t m_nsite;
    size_t m_nelec;

    Metadata() : Archivable("metadata", true, true, true) {

    }
};

/**
 * creates and manages the lifetimes of HDF5 archives
 */
struct Archive {
    const fciqmc_config::Document m_opts;
    const bool m_do_load, m_do_save, m_do_chkpts;
    size_t m_nchkpt = 0ul;
    size_t m_icycle_last_chkpt = 0ul;
    Timer m_timer;
    /**
     * store members and call their load and save methods polymorphically when required
     */
    std::list<Archivable *> m_members;

    explicit Archive(const fciqmc_config::Document &opts);

    ~Archive();

private:

    void save_metadata(hdf5::FileWriter &fw);

    /**
     * check that metadata in loaded file is compatible with options of the current execution
     */
    void verify_metadata(hdf5::FileReader &fr);

    void save(hdf5::FileWriter& fw);


public:

    void load();

    void save();

    void chkpt(const size_t &icycle);

    void add_member(Archivable &item);

    void add_members() {}
    template<typename ...Args>
    void add_members(Archivable &first, Args &... rest) {
        add_member(first);
        add_members(std::forward<Args &>(rest)...);
    }

};

#endif //M7_ARCHIVABLE_H
