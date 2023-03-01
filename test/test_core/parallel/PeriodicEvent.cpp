//
// Created by rja on 01/03/23.
//

#include "test_core/defs.h"
#include "M7_lib/parallel/PeriodicEvent.h"

TEST(PeriodicEvent, CycleAndTimer) {
    PeriodicEvent pe(0, 6, 1);
    uint_t ievent;
    ASSERT_EQ(pe.get_event(0, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(1, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(4, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    Timer::sleep(1.1);
    ASSERT_EQ(pe.get_event(5, ievent), PeriodicEvent::Timed);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(10, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(11, ievent), PeriodicEvent::Cycles);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(12, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 2ul);
    ASSERT_EQ(pe.get_event(17, ievent), PeriodicEvent::Cycles);
    ASSERT_EQ(ievent, 2ul);
    Timer::sleep(1.1);
    ASSERT_EQ(pe.get_event(18, ievent), PeriodicEvent::Timed);
    ASSERT_EQ(ievent, 3ul);
    ASSERT_EQ(pe.get_event(18, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 4ul);
}

TEST(PeriodicEvent, CycleOnly){
    PeriodicEvent pe(0, 6);
    uint_t ievent;
    ASSERT_EQ(pe.get_event(0, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(1, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(4, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    Timer::sleep(1.1);
    ASSERT_EQ(pe.get_event(10, ievent), PeriodicEvent::Cycles);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(12, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(15, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(16, ievent), PeriodicEvent::Cycles);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(16, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 2ul);
}

TEST(PeriodicEvent, TooEarlyCycle){
    PeriodicEvent pe(100, 6);
    uint_t ievent;
    ASSERT_EQ(pe.get_event(100, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(pe.get_event(101, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(pe.get_event(99, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(200, ievent), PeriodicEvent::Cycles);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(200, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 1ul);
}

TEST(PeriodicEvent, TimerOnly) {
    PeriodicEvent pe(0, 0, 1);
    uint_t ievent;
    ASSERT_EQ(pe.get_event(0, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(1, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(4, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 0ul);
    Timer::sleep(1.1);
    ASSERT_EQ(pe.get_event(5, ievent), PeriodicEvent::Timed);
    ASSERT_EQ(ievent, 0ul);
    ASSERT_EQ(pe.get_event(10, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(12, ievent), PeriodicEvent::NotDue);
    Timer::sleep(1.1);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(18, ievent), PeriodicEvent::Timed);
    ASSERT_EQ(ievent, 1ul);
    ASSERT_EQ(pe.get_event(18, ievent), PeriodicEvent::NotDue);
    ASSERT_EQ(ievent, 2ul);
}
