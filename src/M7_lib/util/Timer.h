//
// Created by rja on 02/07/2020.
//

#ifndef M7_TIMER_H
#define M7_TIMER_H

#include <chrono>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>

class Timer {
    std::chrono::duration<double> m_start, m_total;
public:

    Timer();

    bool paused() const;

    void reset();

    void unpause();

    void pause();

    /**
     * @return
     * if paused, the duration between last reset and last pause in seconds.
     * else, the time elapsed since the last reset;
     */
    operator double() const;
};


#endif //M7_TIMER_H
