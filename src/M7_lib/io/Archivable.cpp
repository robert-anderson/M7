//
// Created by Robert J. Anderson on 05/07/2021.
//

#include <M7_lib/hdf5/Group.h>
#include "Archivable.h"

Archivable::Archivable(std::string name, bool load, bool save, bool chkpt) :
        m_name(std::move(name)), m_load(load), m_save(save), m_chkpt(chkpt) {}

Archivable::Archivable(std::string name, const conf::Archivable& opts) :
        Archivable(std::move(name), opts.m_load, opts.m_save, opts.m_chkpt) {}

Archivable::Archivable(std::string name, const conf::Archive& opts) :
        Archivable(std::move(name), opts.do_load(), opts.do_save(), opts.do_chkpts()) {}

void Archivable::load(const hdf5::NodeReader& parent) {
    if (m_load) {
        REQUIRE_TRUE_ALL(parent.child_exists(m_name),
                         log::format("Load failed: \"{}\" does not exist in archive", m_name));
        log::info("loading \"{}\" from archive", m_name);
        load_fn(parent);
    }
}

void Archivable::save(const hdf5::NodeWriter& parent) {
    if (m_save) {
        log::info("saving \"{}\" to archive", m_name);
        save_fn(parent);
    }
}

void Archivable::chkpt(const hdf5::NodeWriter& parent) {
    if (m_chkpt) {
        log::info("saving \"{}\" checkpoint to archive", m_name);
        save_fn(parent);
    }
}

Archive::Archive(const conf::Document& opts) :
        m_opts(opts), m_do_load(m_opts.m_archive.do_load()),
        m_do_save(opts.m_archive.do_save()), m_do_chkpts(opts.m_archive.do_chkpts()) {
    m_timer.unpause();
    if (m_do_load) log::info("reading archive from file \"{}\"", m_opts.m_archive.m_load_path.get());
    else log::info("not reading archive from file");
    if (m_do_save) log::info("saving archive to file \"{}\" upon termination of the solver loop", m_opts.m_archive.m_save_path.get());
    else log::info("not saving archive to file");
    if (m_do_chkpts) {
        log::info("dumping checkpoint archives on file \"{}\"", m_opts.m_archive.m_chkpt_path.get());
        if (m_opts.m_archive.m_period)
            log::info("checkpoints will be made every {} cycles", m_opts.m_archive.m_period.get());
        if (m_opts.m_archive.m_period_mins)
            log::info("checkpoints will be made every {} minutes", m_opts.m_archive.m_period_mins.get());
    } else log::info("not dumping checkpoint archives");
}

Archive::~Archive() {
    save();
}

void Archive::save_metadata(const hdf5::FileWriter& fw) {
    hdf5::GroupWriter gw(fw, "metadata");
    gw.save("timestamp", std::round(std::chrono::steady_clock::now().time_since_epoch().count()));
}

void Archive::verify_metadata(const hdf5::FileReader& /*fr*/) {

}

void Archive::save(const hdf5::FileWriter& fw) {
    save_metadata(fw);
    hdf5::GroupWriter gw(fw, "archive");
    for (auto ptr: m_members) ptr->save(gw);
}

void Archive::load() {
    if (!m_do_load) return;
    hdf5::FileReader fr(m_opts.m_archive.m_load_path);
    verify_metadata(fr);
    hdf5::GroupReader gr(fr, "archive");
    for (auto ptr: m_members) ptr->load(gr);
}

void Archive::save() {
    if (!m_do_save) return;
    hdf5::FileWriter fw(m_opts.m_archive.m_save_path);
    save(fw);
}

void Archive::chkpt(uint_t icycle) {
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

void Archive::add_member(Archivable& item) {
    m_members.push_back(&item);
}
