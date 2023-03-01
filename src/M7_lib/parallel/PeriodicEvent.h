//
// Created by rja on 28/02/23.
//

#ifndef M7_PERIODICEVENT_H
#define M7_PERIODICEVENT_H

#include "M7_lib/conf/Conf.h"
#include "M7_lib/util/Timer.h"
#include "Epoch.h"

/**
 * it is desirable for some processes to be performed periodically, either in terms of wall time, cycle number, or both
 * if both period definitions are in use, the event is triggered by number of seconds or number of periods elapsed since
 * the last event - whichever occurs first
 */
class PeriodicEvent {
    Timer m_timer;
    /**
     * number of cycles in a period
     */
    const uint_t m_ncycle;
    /**
     * number of seconds in a period
     */
    const uint_t m_nsec;
    /**
     * number of events raised so far
     */
    uint_t m_nevent = 0ul;
    /**
     * seconds from beginning of timer to the last event raised, whether triggered by cycles or time
     */
    double m_nsec_last_event = 0.0;
protected:
    /**
     * cycle number of the last event, whether triggered by cycles or time
     */
    uint_t m_icycle_last_event;

public:
    PeriodicEvent(uint_t icycle_begin, uint_t ncycle=0ul, uint_t nsec=0ul);

    enum Reason {NotDue, Timed, Cycles};

    /**
     * @param icycle
     *  cycle number
     * @param ievent
     *  reference to current event index
     * @return
     *  true if event is due
     */
    Reason get_event(uint_t icycle, uint_t& ievent);

    operator bool () const {
        return m_nsec || m_ncycle;
    }

protected:
    virtual Reason due(uint_t icycle);

};

class PeriodicFileSeries : public PeriodicEvent {
    const Epoch& m_epoch;
    const str_t m_path_fmt;

    Reason due(uint_t icycle) override;

public:
    PeriodicFileSeries(const Epoch& epoch, const conf::OptionalFileSeries& series);

    Reason get_file_path(uint_t icycle, str_t& path);

    str_t get_file_path(uint_t icycle);
};


#endif //M7_PERIODICEVENT_H
