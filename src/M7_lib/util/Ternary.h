//
// Created by Robert J. Anderson on 09/10/2020.
//

#ifndef M7_TERNARY_H
#define M7_TERNARY_H

#include <stdexcept>
#include <M7_lib/defs.h>

/*
 * "Ternary" boolean type with a default unspecified state
 */

struct Tern {
    enum State {
        False, True, Neither
    };
private:
    State m_state = Neither;
public:
    Tern(){}
    Tern(bool b){m_state=static_cast<State>(b);}
    Tern(int i){m_state=(i==0||i==1)?static_cast<State>(i):Neither;}
    Tern(const Tern& other){m_state=other.m_state;}
    Tern& operator=(const Tern& other){m_state=other.m_state; return *this;}
    Tern& operator=(bool b){m_state=static_cast<State>(b); return *this;}
    bool operator==(const State& other) const{return m_state==other;}
    bool operator!=(const State& other) const{return m_state!=other;}
    operator bool() const {
        ASSERT(*this!=Neither);
        return m_state;
    }
};

#endif //M7_TERNARY_H
