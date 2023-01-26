//
// Created by Robert J. Anderson on 26/10/2020.
//

#include "Buffer.h"
#include "TableBase.h"
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
    DEBUG_ASSERT_TRUE(m_begin, "this is an unallocated buffer window");
    DEBUG_ASSERT_TRUE(other.m_begin, "can't assign to an unallocated buffer window");
    const auto nbyte_cpy = std::min(other.m_size, m_size);
    /*
     * allocate enough memory for the copied records if needed
     */
    if (m_size < nbyte_cpy) Buffer::Window::resize(nbyte_cpy);
    if (i_can_modify()) std::memcpy(m_begin, other.m_begin, nbyte_cpy);
    m_size = nbyte_cpy;
    m_end = m_begin + m_size;
    m_nrow = m_size / m_row_size;
    m_hwm = m_begin + other.size_in_use();
    return *this;
}

bool Buffer::Window::operator==(const Buffer::Window& other) const {
    if (m_size != other.m_size) return false;
    if (size_in_use() != other.size_in_use()) return false;
    return std::memcmp(m_begin, other.m_begin, size_in_use()) == 0;
}

bool Buffer::Window::allocated() const {
    return m_buffer && m_buffer->size();
}

void Buffer::Window::clear() {
    if (!allocated()) return;
    if (i_can_modify()) std::memset(m_begin, 0, size_in_use());
    m_hwm = m_begin;
}

void Buffer::Window::move(buf_t *begin, uint_t new_size) {
    DEBUG_ASSERT_TRUE(begin, "moving to invalid buffer pointer");
    const auto nbyte_hwm = std::distance(m_begin, m_hwm);
    if (m_begin && i_can_modify()) std::memmove(begin, m_begin, std::min(new_size, m_size));
    m_begin = begin;
    m_size = new_size;
    m_hwm = m_begin + nbyte_hwm;
    m_end = m_begin + m_size;
    m_nrow = m_size / m_row_size;
}

void Buffer::Window::resize(uint_t size, double factor) {
    DEBUG_ASSERT_FALSE(m_buffer==nullptr, "resizing window not associated with buffer");
    m_buffer->resize(size * m_buffer->m_nwindow_max, factor);
}

str_t Buffer::Window::name() const {
    if (!m_buffer) return "";
    return m_buffer->m_name;
}

double Buffer::Window::get_expansion_factor() const {
    return m_buffer->m_expansion_factor;
}

Buffer::Buffer(str_t name, uint_t nwindow_max, bool node_shared) :
        m_name(std::move(name)), m_nwindow_max(nwindow_max), m_node_shared(node_shared) {
    if (!name.empty()) logging::info_("Creating {} buffer \"{}\"",
                                      node_shared ? "node-shared" : "rank-private", name);
    REQUIRE_TRUE(nwindow_max, "A buffer must allow at least one window");
    m_windows.reserve(m_nwindow_max);
}

uint_t Buffer::size() const {
    return m_size;
}

uint_t Buffer::window_size() const {
    return size() / m_nwindow_max;
}

void Buffer::append_window(Buffer::Window *window) {
    REQUIRE_LT(m_windows.size(), m_nwindow_max, "Buffer is over-subscribed");
    if (size()) {
        window->m_begin = m_data + window_size() * m_windows.size();
        window->m_size = window_size();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(uint_t new_size, double factor) {
    if (new_size >= this->size()) {
        // expanding the buffer
        if (factor < 0.0) factor = m_expansion_factor;
    }
    else{
        // shrinking the buffer
        factor = 0.0;
    }
    new_size*= 1.0 + factor;
    // always allocate an integral number of words
    new_size = integer::divceil(new_size, Buffer::c_nbyte_word) * Buffer::c_nbyte_word;
    // quit if the new buffer new_size is the same as the old buffer new_size
    if (new_size == this->size()) return;
    DEBUG_ASSERT_TRUE(new_size, "New size must be non-zero");
    if (!m_name.empty()) {
        logging::info_("Reallocating buffer \"{}\" {} -> {}",
                   m_name, capacity_string(), capacity_string(new_size));
    }

    auto move_windows_fn = [&](buf_t* new_data) {
        auto new_window_size = new_size / m_windows.size();
        for (uint_t iwindow = 0ul; iwindow < m_windows.size(); ++iwindow) {
            // work forwards for shrinking, backwards for enlargement
            auto jwindow = new_size < this->size() ? iwindow : (m_windows.size() - iwindow - 1);
            auto window = m_windows[jwindow];
            DEBUG_ASSERT_FALSE(window->m_buffer == nullptr, "window is not associated with any buffer");
            DEBUG_ASSERT_TRUE(window->m_buffer == this, "window is not associated with this buffer");
            auto new_dbegin = new_data + jwindow * new_window_size;
            window->move(new_dbegin, new_window_size);
        }
    };

    buf_t* tmp_ptr = nullptr;
    /*
     * handle the shared and private cases separately
     */
    if (m_node_shared) {
        SharedArrayBase tmp(new_size, 1);
        tmp_ptr = tmp.m_data;
        move_windows_fn(tmp_ptr);
        m_data_shared = std::move(tmp);
        m_data = m_data_shared.m_data;
    }
    else {
        v_t<buf_t> tmp;
        try {
            tmp.resize(new_size, 0);
        }
        catch (const std::bad_alloc &e) {
            logging::error_("bad allocation");
            ABORT(logging::format("could not allocate sufficient memory to resize buffer \"{}\"", m_name));
        }
        tmp_ptr = tmp.data();
        move_windows_fn(tmp_ptr);
        m_data_priv = std::move(tmp);
        m_data = m_data_priv.data();
        DEBUG_ASSERT_FALSE(m_data_priv.empty(), "new data buffer should be non-empty");
    }
    DEBUG_ASSERT_TRUE(m_data, "new data pointer should be non-null");
    m_size = new_size;
    DEBUG_ASSERT_EQ(m_data, tmp_ptr, "new base ptr was not preserved in move");
    DEBUG_ASSERT_EQ(m_data, m_windows[0]->m_begin, "first window not pointing at begin of buffer");
}

str_t Buffer::capacity_string(uint_t size) const {
    const auto ntable = " x " + std::to_string(m_windows.size()) + " tables";
    if (m_windows.size() == 1)
        return string::memsize(size);
    else
        return string::memsize(size / m_nwindow_max) + ntable + " = " + string::memsize(size);
}

str_t Buffer::capacity_string() const {
    return capacity_string(size());
}