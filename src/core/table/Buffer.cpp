//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"
#include "src/core/util/utils.h"

Buffer::Buffer(std::string name, size_t row_dsize, size_t nrow) :
    m_name(name) {
    resize(row_dsize, nrow);
}

Buffer::~Buffer() {
    //if (!m_name.empty() && dsize()) std::cout << freeing_string() << '\n';
}

Buffer::Buffer(const Buffer &other) : Buffer(other.m_name, other.m_row_dsize, other.m_nrow){}

Buffer::Buffer(Buffer &&other) : Buffer(other.m_name, other.m_row_dsize, other.m_nrow){}

Buffer &Buffer::operator=(const Buffer &other) {
    resize(other.m_row_dsize, other.m_nrow);
    std::copy(other.m_data.begin(), other.m_data.end(), m_data.begin());
    return *this;
}

Buffer &Buffer::operator=(Buffer &&other) {
    m_row_dsize = other.m_row_dsize;
    m_nrow = other.m_nrow;
    //if (!m_name.empty() && dsize()) std::cout << freeing_string() << '\n';
    m_data = std::move(other.m_data);
    return *this;
}

defs::data_t *Buffer::ptr() { return m_data.data(); }

const defs::data_t *Buffer::ptr() const { return m_data.data(); }

size_t Buffer::dsize() const {
    return m_data.size();
}

std::string Buffer::freeing_string() const {
    return "Freeing buffer \"" + m_name + "\" (" + string_utils::memsize(dsize()*defs::nbyte_data)+")";
}

std::string Buffer::capacity_string() const {
    return "Capacity is: " + std::to_string(m_nrow) + " rows of " + std::to_string(m_row_dsize) +
            " datawords each (" + string_utils::memsize(m_row_dsize*m_nrow*defs::nbyte_data) + " in total)";
}

void Buffer::resize(size_t row_dsize, size_t nrow) {
    m_row_dsize = row_dsize;
    m_nrow = nrow;
    if (!m_name.empty() && m_row_dsize && m_nrow) {
        std::cout << "Allocating buffer \"" << m_name <<"\". " << capacity_string()+"\n";
    }
    m_data.resize(m_row_dsize*m_nrow, 0);
    ASSERT(dsize()==m_row_dsize*m_nrow)
}