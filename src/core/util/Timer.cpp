//
// Created by rja on 02/07/2020.
//

#include "Timer.h"

void Timer::unpause() {
    DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when unpausing");
    m_start = std::chrono::steady_clock::now().time_since_epoch();
}

void Timer::pause() {
    DEBUG_ASSERT_FALSE(paused(), "Timer should not be paused when pausing");
    auto now = std::chrono::steady_clock::now().time_since_epoch();
    m_total += now-m_start;
    m_start=std::chrono::duration<double>::max();
}

void Timer::reset() {
    DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when reading value");
    m_total = std::chrono::duration<double>::zero();
}

inline bool Timer::paused() const {
    return m_start==std::chrono::duration<double>::max();
}

Timer::operator double() const {
    DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when reading value");
    return m_total.count();
}