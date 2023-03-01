//
// Created by rja on 28/02/23.
//

#include "PeriodicEvent.h"

PeriodicEvent::PeriodicEvent(uint_t icycle_begin, uint_t ncycle, uint_t nsec) :
        m_ncycle(ncycle), m_nsec(nsec), m_icycle_last_event(icycle_begin){
    m_timer.unpause();
}

PeriodicEvent::Reason PeriodicEvent::due(uint_t icycle) {
    if (icycle<=m_icycle_last_event) return NotDue;
    if (m_timer==m_nsec_last_event) return NotDue;
    if (m_ncycle) {
        // cycles period is in use
        if (m_icycle_last_event + m_ncycle <= icycle) return Cycles;
    }
    if (m_nsec) {
        // wall time period is in use
        // beware: this boolean is not guaranteed to evaluate the same over all MPI ranks
        double time_elapsed = m_timer;
        if (m_nsec_last_event + m_nsec <= time_elapsed) return Timed;
    }
    return NotDue;
}

PeriodicEvent::Reason PeriodicEvent::get_event(uint_t icycle, uint_t& ievent) {
    ievent = m_nevent;
    Reason reason = NotDue;
    auto local_reason = this->due(icycle);
    char is_cycles = (local_reason == Cycles);
    if (mpi::all_land(is_cycles)) reason = Cycles;
    if (reason == NotDue) {
        char is_timed = (local_reason == Timed);
        if (mpi::all_land(is_timed)) reason = Timed;
    }
    if (reason == NotDue) return NotDue;
    // all MPI ranks agree that an event should be raised
    m_nsec_last_event = m_timer;
    m_icycle_last_event = icycle;
    ++m_nevent;
    return reason;
}

PeriodicEvent::Reason PeriodicFileSeries::due(uint_t icycle) {
    if(m_epoch) PeriodicEvent::m_icycle_last_event = m_epoch.icycle_start();
    return PeriodicEvent::due(icycle);
}

PeriodicFileSeries::PeriodicFileSeries(const Epoch& epoch, const conf::OptionalFileSeries& series) :
    PeriodicEvent(~0ul,
        (series.m_enabled && series.m_mode.m_value=="cycle") ? series.m_period.m_value : 0ul,
        (series.m_enabled && series.m_mode.m_value=="minute") ? series.m_period.m_value*60 : 0ul),
    m_epoch(epoch), m_path_fmt(series.m_path_fmt){}

PeriodicEvent::Reason PeriodicFileSeries::get_file_path(uint_t icycle, str_t& path) {
    uint_t ifile;
    auto res = get_event(icycle, ifile);
    path.clear();
    if (res != NotDue) path = logging::format(m_path_fmt, ifile);
    return res;
}

str_t PeriodicFileSeries::get_file_path(uint_t icycle) {
    str_t path;
    get_file_path(icycle, path);
    return path;
}
