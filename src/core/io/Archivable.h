//
// Created by rja on 05/07/2021.
//

#ifndef M7_ARCHIVABLE_H
#define M7_ARCHIVABLE_H

#include <src/core/config/FciqmcConfig.h>
#include <src/core/util/Timer.h>
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

    Archivable(std::string name, bool load, bool save, bool chkpt) :
            m_name(name), m_load(load), m_save(save), m_chkpt(chkpt) {}

    Archivable(std::string name, const fciqmc_config::Io &io_opts) :
            Archivable(name, io_opts.m_load, io_opts.m_save, io_opts.m_chkpt) {}

    /**
     * ctor for object which are not optionally archivable
     * @param name
     *  label of group in HDF5 file
     * @param archive_opts
     *  Global options for the archive
     */
    Archivable(std::string name, const fciqmc_config::Archive &archive_opts) :
            Archivable(name, archive_opts.do_load(), archive_opts.do_save(), archive_opts.do_chkpts()) {}

    virtual ~Archivable() {}

    virtual void load_fn(hdf5::GroupReader &parent) = 0;

    virtual void save_fn(hdf5::GroupWriter &parent) = 0;

    void load(hdf5::GroupReader &parent) {
        if (m_load) {
            REQUIRE_TRUE_ALL(parent.child_exists(m_name),
                             log::format("Load failed: \"{}\" does not exist in archive", m_name));
            log::info("loading \"{}\" from archive", m_name);
            load_fn(parent);
        }
    }

    void save(hdf5::GroupWriter &parent) {
        if (m_save) {
            log::info("saving \"{}\" to archive", m_name);
            save_fn(parent);
        }
    }

    void chkpt(hdf5::GroupWriter &parent) {
        if (m_chkpt) {
            log::info("saving \"{}\" checkpoint to archive", m_name);
            save_fn(parent);
        }
    }
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
    const fciqmc_config::Document &m_opts;
    const bool m_do_load, m_do_save, m_do_chkpts;
    size_t m_nchkpt = 0ul;
    size_t m_icycle_last_chkpt = 0ul;
    Timer m_timer;
    std::list<Archivable *> m_contents;


    Archive(const fciqmc_config::Document &opts) :
            m_opts(opts), m_do_load(m_opts.m_archive.do_load()),
            m_do_save(opts.m_archive.do_save()), m_do_chkpts(opts.m_archive.do_chkpts()) {
        m_timer.unpause();
        if (m_do_load) log::info("reading archive from file \"{}\"", m_opts.m_archive.m_load_path.get());
        else log::info("not reading archive from file");
        if (m_do_save) log::info("saving archive to file \"{}\"", m_opts.m_archive.m_load_path.get());
        else log::info("not saving archive to file");
        if (m_do_chkpts) {
            log::info("dumping checkpoint archives on file \"{}\"", m_opts.m_archive.m_chkpt_path.get());
            if (m_opts.m_archive.m_period)
                log::info("checkpoints will be made every {} cycles", m_opts.m_archive.m_period.get());
            if (m_opts.m_archive.m_period_mins)
                log::info("checkpoints will be made every {} minutes", m_opts.m_archive.m_period_mins.get());
        } else log::info("not dumping checkpoint archives");
    }

private:

    void save_metadata(hdf5::FileWriter &fw) {
        hdf5::GroupWriter gw("metadata", fw);
        gw.save("timestamp", std::chrono::steady_clock::now().time_since_epoch().count());
    }

    /**
     * check that metadata in loaded file is compatible with options of the current execution
     */
    void verify_metadata(hdf5::FileReader &fr) {

    }

    void save(hdf5::FileWriter fw) {
        save_metadata(fw);
        hdf5::GroupWriter gw("archive", fw);
        for (auto ptr: m_contents) ptr->save(gw);
    }


public:

    void load() {
        if (m_opts.m_archive.m_load_path.get().empty()) return;
        hdf5::FileReader fr(m_opts.m_archive.m_load_path);
        verify_metadata(fr);
        hdf5::GroupReader gr("archive", fr);
        for (auto ptr: m_contents) ptr->load(gr);
    }

    void save() {
        if (m_opts.m_archive.m_save_path.get().empty()) return;
        hdf5::FileWriter fw(m_opts.m_archive.m_save_path);
        save(fw);
    }

    void chkpt(const size_t &icycle) {
        if (m_do_chkpts) return;
        bool output = false;
        output = m_opts.m_archive.m_period_mins && m_timer / 60 >= m_opts.m_archive.m_period_mins;
        if (output) m_timer.reset();
        else if (m_opts.m_archive.m_period)
            output = (icycle - m_icycle_last_chkpt) == m_opts.m_archive.m_period;
        // never output chkpt twice on the same cycle
        if (output && (icycle == m_icycle_last_chkpt)) output = false;
        if (!output) return;

        hdf5::FileWriter fw(log::format(m_opts.m_archive.m_chkpt_path, m_nchkpt));
        save(fw);
    }

    void add_content(Archivable *item) {
        m_contents.push_back(item);
    }

};

#endif //M7_ARCHIVABLE_H
