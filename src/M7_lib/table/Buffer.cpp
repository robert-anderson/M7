//
// Created by Robert J. Anderson on 26/10/2020.
//

#include "Buffer.h"
#include <M7_lib/io/Logging.h>
#include <M7_lib/parallel/MPIAssert.h>
#include "M7_lib/util/String.h"
#include <M7_lib/util/Integer.h>


Buffer::Window::Window(Buffer *buffer, uint_t row_size): Window(row_size) {
    ASSERT(buffer);
    buffer->append_window(this);
    ASSERT(m_buffer == buffer);
    ASSERT(m_buffer->m_windows.back() == this);
}

Buffer::Window &Buffer::Window::operator=(const Buffer::Window &other) {
    DEBUG_ASSERT_EQ(other.m_row_size, m_row_size, "can't assign to incompatible window");
    DEBUG_ASSERT_FALSE(m_begin==nullptr, "this is an unallocated buffer window");
    DEBUG_ASSERT_FALSE(other.m_begin==nullptr, "can't assign to an unallocated buffer window");
    auto nbyte = std::min(other.m_size, m_size);
    std::memcpy(m_begin, other.m_begin, nbyte);
    m_size = nbyte;
    m_nrow = m_size/m_row_size;
    return *this;
}

bool Buffer::Window::allocated() const {
    return m_buffer && m_buffer->size();
}

void Buffer::Window::move(defs::buf_t *begin, uint_t new_size) {
    DEBUG_ASSERT_FALSE(begin==nullptr, "moving to invalid buffer pointer");
    if (m_begin) std::memmove(begin, m_begin, m_size);
    m_begin = begin;
    m_size = new_size;
    m_nrow = m_size/m_row_size;
}

void Buffer::Window::resize(uint_t size, double factor) {
    DEBUG_ASSERT_FALSE(m_buffer==nullptr, "resizing window not associated with buffer");
    m_buffer->resize(size * m_buffer->m_nwindow_max, factor);
}

std::string Buffer::Window::name() const {
    if (!m_buffer) return "";
    return m_buffer->m_name;
}

double Buffer::Window::get_expansion_factor() const {
    return m_buffer->m_expansion_factor;
}


Buffer::Buffer(std::string name, uint_t nwindow_max) :
        m_name(std::move(name)), m_nwindow_max(nwindow_max) {
    if (!name.empty()) log::info_("Creating \"{}\" buffer", name);
    REQUIRE_TRUE(nwindow_max, "A buffer must allow at least one window");
    m_windows.reserve(m_nwindow_max);
}

uint_t Buffer::size() const {
    return m_data.size();
}

uint_t Buffer::window_size() const {
    return size() / m_nwindow_max;
}

void Buffer::append_window(Buffer::Window *window) {
    REQUIRE_LT(m_windows.size(), m_nwindow_max, "Buffer is over-subscribed");
    if (size()) {
        window->m_begin = m_data.data() + window_size() * m_windows.size();
        window->m_size = window_size();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(uint_t size, double factor) {
    if (size>=this->size()) {
        // expanding the buffer
        if (factor < 0.0) factor = m_expansion_factor;
    }
    else{
        // shrinking the buffer
        factor = 0.0;
    }
    size*=1.0+factor;
    // always allocate an integral number of words
    size = integer::divceil(size, Buffer::c_nbyte_word)*Buffer::c_nbyte_word;
    DEBUG_ASSERT_TRUE(size, "New size must be non-zero");
    if (!m_name.empty()) {
        log::info_("Reallocating buffer \"{}\" {} -> {}",
                   m_name, capacity_string(), capacity_string(size));
    }
    std::vector<defs::buf_t> tmp;
    try {
        tmp.resize(size, 0);
    }
    catch (const std::bad_alloc& e){
        ABORT(log::format("could not allocate sufficient memory to resize buffer \"{}\"", m_name));
    }
    auto new_window_size = tmp.size() / m_nwindow_max;

    for (uint_t iwindow = 0ul; iwindow < m_windows.size(); ++iwindow) {
        // work forwards for shrinking, backwards for enlargement
        auto jwindow = size < this->size() ? iwindow : m_windows.size() - iwindow - 1;
        auto window = m_windows[jwindow];
        DEBUG_ASSERT_FALSE(window->m_buffer==nullptr, "window is not associated with any buffer");
        DEBUG_ASSERT_TRUE(window->m_buffer==this, "window is not associated with this buffer");
        auto new_dbegin = tmp.data() + jwindow * new_window_size;
        window->move(new_dbegin, new_window_size);
    }
    auto tmp_ptr = tmp.data();
    DEBUG_ONLY(tmp_ptr);
    m_data = std::move(tmp);
    DEBUG_ASSERT_EQ(m_data.data(), tmp_ptr, "new base ptr was not preserved in move");
    DEBUG_ASSERT_EQ(m_data.size(), this->size(), "new size is incorrect");
    DEBUG_ASSERT_EQ(m_data.data(), m_windows[0]->m_begin, "first window not pointing at begin of buffer");
}

std::string Buffer::capacity_string(uint_t size) const {
    const auto ntable = " x " + std::to_string(m_windows.size()) + " tables";
    if (m_windows.size() == 1)
        return string::memsize(size);
    else
        return string::memsize(size / m_nwindow_max) + ntable + " = " + string::memsize(size);
}

std::string Buffer::capacity_string() const {
    return capacity_string(size());
}
