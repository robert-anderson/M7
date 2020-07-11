//
// Created by rja on 02/07/2020.
//

#ifndef M7_TIMER_H
#define M7_TIMER_H

#include <chrono>
#include "defs.h"

class Timer {
    std::chrono::duration<double> m_start, m_lap, m_total;
public:

    Timer():m_start{std::chrono::duration<double>::max()} {}

    void reset();

    void unpause();

    void pause();

    double total();
    double lap();

};


#endif //M7_TIMER_H
