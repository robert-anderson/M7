//
// Created by rja on 02/07/2020.
//

#ifndef M7_TIMER_H
#define M7_TIMER_H

#include <chrono>
#include "src/defs.h"
#include "src/core/parallel/MPIAssert.h"

class Timer {
    std::chrono::duration<double> m_start, m_total;
public:

    Timer():m_start{std::chrono::duration<double>::max()} {}

    bool paused() const{
        return m_start==std::chrono::duration<double>::max();
    }

    void reset() {
        DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when reading value");
        m_total = std::chrono::duration<double>::zero();
    }

    void unpause() {
        DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when unpausing");
        m_start = std::chrono::steady_clock::now().time_since_epoch();
    }

    void pause() {
        DEBUG_ASSERT_FALSE(paused(), "Timer should not be paused when pausing");
        auto now = std::chrono::steady_clock::now().time_since_epoch();
        m_total += now-m_start;
        m_start=std::chrono::duration<double>::max();
    }

    /**
     * @return
     * if paused, the duration between last reset and last pause in seconds.
     * else, the time elapsed since the last reset;
     */
    operator double() const {
        if (paused()) return m_total.count();
        else return (std::chrono::steady_clock::now().time_since_epoch() - m_start).count();
    }
};


#endif //M7_TIMER_H
