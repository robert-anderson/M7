//
// Created by rja on 02/07/2020.
//

#include "Timer.h"

void Timer::unpause() {
    ASSERT(m_start==std::chrono::duration<double>::max())
    m_start = std::chrono::steady_clock::now().time_since_epoch();
}

std::chrono::duration<double> Timer::pause() {
    ASSERT(m_start!=std::chrono::duration<double>::max())
    auto now = std::chrono::steady_clock::now().time_since_epoch();
    auto diff = now-m_start;
    m_elapsed += diff;
    m_start=std::chrono::duration<double>::max();
    return diff;
}

double Timer::elapsed() {
    return m_elapsed.count();
}
