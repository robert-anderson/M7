//
// Created by Robert J. Anderson on 02/07/2020.
//

#include "Timer.h"
#include <M7_lib/parallel/MPIAssert.h>


void Timer::sleep_(double nsec) {
    Timer t;
    t.unpause();
    while (t<nsec){}
}

void Timer::sleep(double nsec) {
    mpi::barrier();
    sleep_(nsec);
    mpi::barrier();
}

Timer::Timer() :m_start{std::chrono::duration<double>::max()} {}

bool Timer::paused() const {
    return m_start==std::chrono::duration<double>::max();
}

void Timer::reset() {
    DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when reading value");
    m_total = std::chrono::duration<double>::zero();
}

void Timer::unpause() {
    DEBUG_ASSERT_TRUE(paused(), "Timer should be paused when unpausing");
    m_start = std::chrono::steady_clock::now().time_since_epoch();
}

void Timer::pause() {
    DEBUG_ASSERT_FALSE(paused(), "Timer should not be paused when pausing");
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    m_total += now-m_start;
    m_start=std::chrono::duration<double>::max();
}

Timer::operator double() const {
    if (paused()) return m_total.count();
    else {
        const auto now = std::chrono::steady_clock::now().time_since_epoch();
        const std::chrono::duration<double> diff = now-m_start;
        return diff.count();
    }
}