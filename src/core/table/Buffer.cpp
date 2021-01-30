//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"
#include "src/core/io/Logging.h"
#include "src/core/parallel/MPIAssert.h"

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

void Buffer::Window::make_room(size_t nrow) {
    ASSERT(m_buffer)
    m_buffer->make_room(nrow * m_buffer->m_nwindow_max);
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
    if (!name.empty()) log::info_("Creating \"{}\" buffer", name);
    MPI_REQUIRE(nwindow_max,"A buffer must allow at least one window");
    MPI_REQUIRE(row_dsize, "The row must consist of a non-zero number datawords");
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


void Buffer::append_window(Buffer::Window *window) {
    MPI_REQUIRE_ALL(m_windows.size() < m_nwindow_max, "Buffer is over-subscribed");
    if (dsize()) {
        window->m_dbegin = m_data.data()+window_dsize() * m_windows.size();
        window->m_dend = window->m_dbegin + window_dsize();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(size_t nrow) {
    if (!m_name.empty()) {
        log::info_("Reallocating buffer \"{}\" {} -> {}",
        m_name, capacity_string(), capacity_string(nrow*m_row_dsize));
    }
    std::vector<defs::data_t> tmp(nrow*m_row_dsize, 0ul);
    auto new_window_dsize = tmp.size() / m_nwindow_max;

    for (size_t iwindow = 0ul; iwindow < m_windows.size(); ++iwindow) {
        // work backwards for enlargement
        auto jwindow = m_windows.size() - iwindow - 1;
        auto window = m_windows[jwindow];
        ASSERT(window->m_buffer)
        ASSERT(window->m_buffer==this)
        auto new_dbegin = tmp.data() + jwindow * new_window_dsize;
        window->move(new_dbegin, new_dbegin + new_window_dsize);
    }
    m_data = std::move(tmp);
    ASSERT(m_data.data() ==m_windows[0]->dbegin());
}

void Buffer::make_room(size_t nrow) {
    if (dsize()<m_row_dsize*nrow) resize(nrow);
}

void Buffer::expand(size_t delta_nrow, double expansion_factor) {
    MPI_REQUIRE(expansion_factor>0.0, "invalid expansion factor");
    resize((nrow() + delta_nrow)*(1+expansion_factor));
}

void Buffer::expand(size_t delta_nrow) {
    expand(delta_nrow, m_expansion_factor);
}

std::string Buffer::capacity_string(size_t dsize) const {
    const auto ntable = "x "+std::to_string(m_windows.size())+" tables";
    const auto per_table = std::to_string(dsize/(m_row_dsize*m_nwindow_max));
    const auto total = std::to_string(dsize/m_row_dsize);
    if (m_windows.size()==1)
        return "[" + per_table +" rows (" + string_utils::memsize(dsize*defs::nbyte_data) + ")]";
    else
        return "[" + per_table +" rows "+ntable+" = " + total + " total (" + string_utils::memsize(dsize*defs::nbyte_data) + ")]";
}

std::string Buffer::capacity_string() const {
    return capacity_string(dsize());
}