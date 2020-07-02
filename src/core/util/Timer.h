//
// Created by rja on 02/07/2020.
//

#ifndef M7_TIMER_H
#define M7_TIMER_H

#include <chrono>
#include "defs.h"

class Timer {
    std::chrono::duration<double> m_start, m_elapsed;
public:

    Timer():m_start{std::chrono::duration<double>::max()} {}

    void unpause();

    std::chrono::duration<double> pause();

    double elapsed();

};


#endif //M7_TIMER_H
