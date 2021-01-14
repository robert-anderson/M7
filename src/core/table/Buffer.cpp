//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"


Buffer::Window::Window(Buffer *buffer) {
    ASSERT(buffer);
    buffer->append_window(this);
    ASSERT(m_buffer==buffer);
    ASSERT(m_buffer->m_windows.back()==this);
}

size_t Buffer::Window::dsize() const {
    return std::distance(m_dbegin, m_dend);
}

size_t Buffer::Window::size() const {
    return defs::nbyte_data * dsize();
}

void Buffer::Window::move(defs::data_t *dbegin, defs::data_t *dend) {
    if (m_dbegin) std::memmove(dbegin, m_dbegin, size());
    m_dbegin = dbegin;
    m_dend = dend;
}

defs::data_t *Buffer::Window::dbegin() {
    return m_dbegin;
}

const defs::data_t *Buffer::Window::dbegin() const {
    return m_dbegin;
}

void Buffer::Window::resize(size_t dsize) {
    ASSERT(m_buffer)
    m_buffer->resize(dsize * m_buffer->m_nwindow_max);
}

void Buffer::Window::expand() {
    ASSERT(m_buffer)
    m_buffer->expand();
}

void Buffer::Window::expand(size_t delta_dsize) {
    ASSERT(m_buffer)
    m_buffer->expand(delta_dsize * m_buffer->m_nwindow_max);
}

Buffer::Buffer(std::string name, size_t nwindow_max) :
        m_name(std::move(name)), m_nwindow_max(nwindow_max) {
    if (!nwindow_max) mpi::stop_all("A buffer must allow at least one window");
    m_windows.reserve(m_nwindow_max);
}

Buffer::Buffer(std::string name, size_t nwindow_max, size_t dsize) :
        Buffer(name, nwindow_max) {
    m_data.resize(dsize, 0);
}

size_t Buffer::dsize() const {
    return m_data.size();
}

size_t Buffer::window_dsize() const {
    return dsize() / m_nwindow_max;
}

defs::data_t *Buffer::dbegin() {
    return m_data.data();
}

const defs::data_t *Buffer::dbegin() const {
    return m_data.data();
}

defs::data_t *Buffer::dbegin(const size_t &iwindow) {
    return m_data.data() + window_dsize() * iwindow;
}

const defs::data_t *Buffer::dbegin(const size_t &iwindow) const {
    return m_data.data() + window_dsize() * iwindow;
}

void Buffer::append_window(Buffer::Window *window) {
    if (m_windows.size() >= m_nwindow_max) mpi::stop_all("Buffer is over-subscribed");
    if (dbegin()) {
        window->m_dbegin = dbegin(m_windows.size());
        window->m_dend = window->m_dbegin + window_dsize();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(size_t dsize) {
    if (!m_name.empty()) {
        std::cout << "Reallocating buffer \"" << m_name << "\". " << capacity_string()+"\n";
    }
    std::vector<defs::data_t> tmp(dsize, 0ul);
    auto new_window_dsize = dsize / m_nwindow_max;
    for (size_t iwindow = 0ul; iwindow < m_nwindow_max; ++iwindow) {
        // work backwards for enlargement
        auto window = m_windows[m_windows.size() - iwindow - 1];
        ASSERT(window->m_buffer)
        ASSERT(window->m_buffer==this)
        auto new_dbegin = tmp.data() + iwindow * new_window_dsize;
        window->move(new_dbegin, new_dbegin + new_window_dsize);
    }
    m_data = std::move(tmp);
}

void Buffer::expand() {
    expand(size_t(dsize()*m_expansion_factor));
}

void Buffer::expand(size_t delta_dsize) {
    resize(dsize() + delta_dsize);
}

std::string Buffer::capacity_string() const {
    return "Capacity is: " + std::to_string(dsize())+" datawords (" + string_utils::memsize(dsize()*defs::nbyte_data) + ")";
}
