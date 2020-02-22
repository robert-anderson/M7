//
// Created by Robert John Anderson on 2020-02-04.
//

#include <assert.h>
#include <cstring>
#include <iostream>
#include <x86intrin.h>
#include "BitfieldNew.h"
#include "../defs.h"
#include "BitfieldHasher.h"

BitfieldNew::BitfieldNew(const size_t &nbit):
m_data_internal(ndataword(nbit)), m_nbit(nbit), m_ndataword(ndataword(nbit)){
    m_data = m_data_internal.data();
}

BitfieldNew::BitfieldNew(const size_t &nbit, defs::data_t *data_external):
m_nbit(nbit), m_ndataword(ndataword(nbit)), m_data_external(data_external){
    m_data = m_data_external;
}

void BitfieldNew::set(const size_t &i) {
    assert(i<m_nbit);
    set_bit(m_data[i/(8*sizeof(defs::data_t))], i%(8*sizeof(defs::data_t)));
}

void BitfieldNew::set(const defs::inds &inds) {
    for (auto i: inds) set(i);
}

void BitfieldNew::clr(const size_t &i) {
    assert(i<m_nbit);
    clr_bit(m_data[i/(8*sizeof(defs::data_t))], i%(8*sizeof(defs::data_t)));
}

bool BitfieldNew::get(const size_t &i) const {
    assert(i<m_nbit);
    return get_bit(m_data[i/(8*sizeof(defs::data_t))], i%(8*sizeof(defs::data_t)));
}

defs::data_t BitfieldNew::get_dataword(const size_t &i) const {
    assert(i<m_ndataword);
    return m_data[i];
}

void BitfieldNew::zero(){
    memset((void*)m_data, 0, m_ndataword*sizeof(defs::data_t));
}

std::string BitfieldNew::to_string(size_t padding) const {
    std::string out{};
    for (size_t i =0ul; i < m_nbit; ++i) {
        out+=get(i)?"1":"0";
    }
    out.insert(out.begin(), padding, ' ');
    return out;
}

void BitfieldNew::print() const {
    std::cout << to_string() << std::endl;
}


size_t BitfieldNew::nsetbits() const {
    size_t n = 0ul;
    for (size_t i=0ul; i < m_ndataword; ++i) {
        if (sizeof(defs::data_t)==8) {
            n += _popcnt64(m_data[i]);
        }
        else if (sizeof(defs::data_t)==4) {
            n += _popcnt32(m_data[i]);
        }
    }
    return n;
}

size_t BitfieldNew::nsetbits_common(const BitfieldNew &other) const {
    size_t n = 0ul;
    for (size_t i=0ul; i < m_ndataword; ++i) {
        if (sizeof(defs::data_t)==8) {
            n += _popcnt64(m_data[i]&other.m_data[i]);
        }
        else if (sizeof(defs::data_t)==4) {
            n += _popcnt32(m_data[i]&other.m_data[i]);
        }
    }
    return n;
}


size_t BitfieldNew::nsetbits_cleared(const BitfieldNew &other) const {
    size_t n = 0ul;
    for (size_t i=0ul; i < m_ndataword; ++i) {
        if (sizeof(defs::data_t)==8) {
            n += _popcnt64(m_data[i]&~other.m_data[i]);
        }
        else if (sizeof(defs::data_t)==4) {
            n += _popcnt32(m_data[i]&~other.m_data[i]);
        }
    }
    return n;
}

size_t BitfieldNew::hash(const size_t modular_divisor) const {
    return BitfieldHasher()(*this)%modular_divisor;
}

bool BitfieldNew::is_zero() const {
    return nsetbits()==0;
}

bool BitfieldNew::is_null() const {
    return m_nbit==0;
}
