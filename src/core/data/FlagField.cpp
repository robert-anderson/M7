//
// Created by rja on 01/11/2020.
//

#include "FlagField.h"

Flag::Flag(FlagField *field, size_t nbit) :
        m_nbit(nbit), m_offset(field->add_flag(this)){}
