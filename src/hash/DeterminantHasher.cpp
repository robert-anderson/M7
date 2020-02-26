//
// Created by Robert John Anderson on 2020-02-13.
//

#include "DeterminantHasher.h"
#include "BitfieldHasher.h"

size_t DeterminantHasher::operator()(const Determinant &key) const {
    return BitfieldHasher()(key.m_bitfields[0])^BitfieldHasher()(key.m_bitfields[1]);
}