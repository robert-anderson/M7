//
// Created by rja on 02/07/2020.
//

#include "Timer.h"

void Timer::unpause() {
    MPI_ASSERT(paused(), "Timer should be paused when unpausing");
    m_start = std::chrono::steady_clock::now().time_since_epoch();
}

void Timer::pause() {
    MPI_ASSERT(!paused(), "Timer should not be paused when pausing");
    auto now = std::chrono::steady_clock::now().time_since_epoch();
    m_total += now-m_start;
    m_start=std::chrono::duration<double>::max();
}

void Timer::reset() {
    MPI_ASSERT(paused(), "Timer should be paused when reading value");
    m_total = std::chrono::duration<double>::zero();
}

inline bool Timer::paused() const {
    return m_start==std::chrono::duration<double>::max();
}

Timer::operator double() const {
    MPI_ASSERT(paused(), "Timer should be paused when reading value");
    return m_total.count();
}