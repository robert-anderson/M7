//
// Created by rja on 02/07/2020.
//

#include "Timer.h"

void Timer::unpause() {
    ASSERT(m_start==std::chrono::duration<double>::max())
    m_start = std::chrono::steady_clock::now().time_since_epoch();
}

void Timer::pause() {
    ASSERT(m_start!=std::chrono::duration<double>::max())
    auto now = std::chrono::steady_clock::now().time_since_epoch();
    m_lap = now-m_start;
    m_total += m_lap;
    m_start=std::chrono::duration<double>::max();
}

double Timer::total() {
    return m_total.count();
}

double Timer::lap() {
    return m_lap.count();
}

void Timer::reset() {
    ASSERT(m_start==std::chrono::duration<double>::max())
    m_lap = std::chrono::duration<double>::zero();
    m_total = std::chrono::duration<double>::zero();
}
