//
// Created by rja on 02/10/2020.
//

#include "Buffer.h"

Buffer::Buffer() : m_data() {}

Buffer::Buffer(size_t dsize) : Buffer() {
    resize(dsize);
}

Buffer::Buffer(size_t ncacheline, size_t nrow) : Buffer() {
    resize(ncacheline * defs::ncacheline_data * nrow);
}

defs::data_t *Buffer::ptr() { return m_data.data(); }

const defs::data_t *Buffer::ptr() const { return m_data.data(); }

size_t Buffer::size() const {
    return m_data.size();
}

void Buffer::resize(size_t dsize) {
    ASSERT(dsize);
    m_data.resize(dsize);
}

void Buffer::resize(size_t ncacheline, size_t nrow) {
    ASSERT(ncacheline);
    resize((nrow * ncacheline) * defs::ncacheline_data);
}
