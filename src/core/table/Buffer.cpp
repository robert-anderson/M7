//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"
#include "src/core/io/Logging.h"

Buffer::Window::Window(Buffer *buffer) {
    ASSERT(buffer);
    buffer->append_window(this);
    ASSERT(m_buffer==buffer);
    ASSERT(m_buffer->m_windows.back()==this);
}

bool Buffer::Window::allocated() const {
    return m_buffer && m_buffer->dsize();
}

size_t Buffer::Window::dsize() const {
    if (!m_dbegin) return 0;
    ASSERT(m_dend);
    ASSERT(m_buffer)
    ASSERT(m_buffer->window_dsize()==(size_t)std::distance(m_dbegin, m_dend));
    return std::distance(m_dbegin, m_dend);
}

size_t Buffer::Window::size() const {
    return defs::nbyte_data * dsize();
}

size_t Buffer::Window::nrow() const {
    return m_buffer->window_nrow();
}

void Buffer::Window::move(defs::data_t *dbegin, defs::data_t *dend) {
    if (m_dbegin) std::memmove(dbegin, m_dbegin, size());
    m_dbegin = dbegin;
    m_dend = dend;
}

defs::data_t *Buffer::Window::dbegin() {
    ASSERT(m_buffer);
    ASSERT(allocated());
    return m_dbegin;
}

const defs::data_t *Buffer::Window::dbegin() const {
    return m_dbegin;
}

void Buffer::Window::resize(size_t nrow) {
    ASSERT(m_buffer)
    m_buffer->resize(nrow * m_buffer->m_nwindow_max);
}

void Buffer::Window::expand(size_t delta_nrow) {
    ASSERT(m_buffer)
    m_buffer->expand(delta_nrow * m_buffer->m_nwindow_max);
}

void Buffer::Window::expand(size_t delta_nrow, double expansion_factor) {
    ASSERT(m_buffer)
    m_buffer->expand(delta_nrow * m_buffer->m_nwindow_max, expansion_factor);
}

double Buffer::Window::expansion_factor() const {
    return m_buffer->m_expansion_factor;
}

Buffer::Buffer(std::string name, size_t nwindow_max, size_t row_dsize) :
        m_name(std::move(name)), m_nwindow_max(nwindow_max), m_row_dsize(row_dsize) {
    if (!nwindow_max) mpi::stop_all("A buffer must allow at least one window");
    if (!row_dsize) mpi::stop_all("The row must consist of a non-zero number datawords");
    m_windows.reserve(m_nwindow_max);
}

size_t Buffer::dsize() const {
    return m_data.size();
}

size_t Buffer::nrow() const {
    return dsize()/m_row_dsize;
}

size_t Buffer::window_dsize() const {
    return dsize() / m_nwindow_max;
}

size_t Buffer::window_nrow() const {
    return nrow() / m_nwindow_max;
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
    if (dsize()) {
        window->m_dbegin = dbegin(m_windows.size());
        window->m_dend = window->m_dbegin + window_dsize();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(size_t nrow) {
    if (!m_name.empty()) {
        log::info("Reallocating buffer \"" + m_name + "\". " + capacity_string());
    }
    std::vector<defs::data_t> tmp(nrow*m_row_dsize, 0ul);
    auto new_window_dsize = tmp.size() / m_nwindow_max;
    for (size_t iwindow = 0ul; iwindow < m_nwindow_max; ++iwindow) {
        // work backwards for enlargement
        auto jwindow = m_windows.size() - iwindow - 1;
        auto window = m_windows[jwindow];
        ASSERT(window->m_buffer)
        ASSERT(window->m_buffer==this)
        auto new_dbegin = tmp.data() + jwindow * new_window_dsize;
        window->move(new_dbegin, new_dbegin + new_window_dsize);
    }
    m_data = std::move(tmp);
    if (!m_name.empty()) {
        log::info("New " + capacity_string());
    }
    ASSERT(dbegin()==m_windows[0]->dbegin());
}

void Buffer::expand(size_t delta_nrow, double expansion_factor) {
    if (expansion_factor<0.0) mpi::stop_all("invalid expansion factor");
    resize((nrow() + delta_nrow)*(1+expansion_factor));
}

void Buffer::expand(size_t delta_nrow) {
    expand(delta_nrow, m_expansion_factor);
}

std::string Buffer::capacity_string() const {
    const auto ntable = "x "+std::to_string(m_windows.size())+" tables";
    const auto per_table = std::to_string(dsize()/(m_row_dsize*m_windows.size()));
    const auto total = std::to_string(dsize()/m_row_dsize);
    if (m_windows.size()==1)
        return "Capacity is: " + per_table +" rows (" + string_utils::memsize(dsize()*defs::nbyte_data) + ")";
    else
        return "Capacity is: " + per_table +" rows "+ntable+" = " + total + " total (" + string_utils::memsize(dsize()*defs::nbyte_data) + ")";
}
