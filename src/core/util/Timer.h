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

    bool paused() const;

    void reset();

    void unpause();

    void pause();

    operator double() const;


};


#endif //M7_TIMER_H
