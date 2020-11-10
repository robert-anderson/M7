//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"

Buffer::Buffer() : m_data() {}

Buffer::Buffer(size_t dsize) : Buffer() {
    resize(dsize);
}

Buffer::Buffer(size_t ndword, size_t nrow) : Buffer() {
    resize(ndword, nrow);
}

defs::data_t *Buffer::ptr() { return m_data.data(); }

const defs::data_t *Buffer::ptr() const { return m_data.data(); }

size_t Buffer::dsize() const {
    return m_data.size();
}

void Buffer::resize(size_t dsize) {
    ASSERT(dsize);
    m_data.resize(dsize);
}

void Buffer::resize(size_t ndword, size_t nrow) {
    ASSERT(ndword);
    resize(nrow * ndword);
}