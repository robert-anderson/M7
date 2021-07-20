//
// Created by RJA on 26/10/2020.
//

#include "Buffer.h"
#include "src/core/io/Logging.h"
#include "src/core/parallel/MPIAssert.h"

Buffer::Window::Window(Buffer *buffer) {
    ASSERT(buffer);
    buffer->append_window(this);
    ASSERT(m_buffer == buffer);
    ASSERT(m_buffer->m_windows.back() == this);
}

bool Buffer::Window::allocated() const {
    return m_buffer && m_buffer->size();
}

size_t Buffer::Window::size() const {
    if (!m_begin) return 0;
    DEBUG_ASSERT_TRUE(m_end, "end pointer is not set");
    DEBUG_ASSERT_TRUE(m_buffer, "buffer is not allocated")
    DEBUG_ASSERT_EQ(m_buffer->window_size(), static_cast<size_t>(std::distance(m_begin, m_end)),
                    "the begin and end pointers of the buffer window are not compatible with its size");
    return std::distance(m_begin, m_end);
}

void Buffer::Window::move(defs::buf_t *begin, defs::buf_t *end) {
    if (m_begin) std::memmove(begin, m_begin, size());
    m_begin = begin;
    m_end = end;
}

void Buffer::Window::resize(size_t size) {
    ASSERT(m_buffer)
    m_buffer->resize(size * m_buffer->m_nwindow_max);
}

std::string Buffer::Window::name() const {
    return m_buffer->m_name;
}

double Buffer::Window::get_expansion_factor() const {
    return m_buffer->m_expansion_factor;
}

Buffer::Buffer(std::string name, size_t nwindow_max) :
        m_name(std::move(name)), m_nwindow_max(nwindow_max) {
    if (!name.empty()) log::info_("Creating \"{}\" buffer", name);
    REQUIRE_TRUE(nwindow_max, "A buffer must allow at least one window");
    m_windows.reserve(m_nwindow_max);
}

size_t Buffer::size() const {
    return m_data.size();
}

size_t Buffer::window_size() const {
    return size() / m_nwindow_max;
}

void Buffer::append_window(Buffer::Window *window) {
    REQUIRE_LT(m_windows.size(), m_nwindow_max, "Buffer is over-subscribed");
    if (size()) {
        window->m_begin = m_data.data() + window_size() * m_windows.size();
        window->m_end = window->m_begin + window_size();
    }
    window->m_buffer = this;
    m_windows.push_back(window);
}

void Buffer::resize(size_t size, double factor) {
    if (size>this->size()) {
        // expanding the buffer
        if (factor < 0.0) factor = m_expansion_factor;
    }
    else{
        // shrinking the buffer
        factor = 0.0;
    }
    size*=1.0+factor;
    // always allocate an integral number of words
    size = integer_utils::divceil(size, defs::nbyte_word)*defs::nbyte_word;
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

    for (size_t iwindow = 0ul; iwindow < m_windows.size(); ++iwindow) {
        // work forwards for shrinking, backwards for enlargement
        auto jwindow = size < this->size() ? iwindow : m_windows.size() - iwindow - 1;
        auto window = m_windows[jwindow];
        ASSERT(window->m_buffer)
        ASSERT(window->m_buffer == this)
        auto new_dbegin = tmp.data() + jwindow * new_window_size;
        window->move(new_dbegin, new_dbegin + new_window_size);
    }
    m_data = std::move(tmp);
    ASSERT(m_data.size() == this->size());
    ASSERT(m_data.data() == m_windows[0]->m_begin);
}

std::string Buffer::capacity_string(size_t size) const {
    const auto ntable = " x " + std::to_string(m_windows.size()) + " tables";
    if (m_windows.size() == 1)
        return "[ " + string_utils::memsize(size) + " ]";
    else
        return "[" + string_utils::memsize(size / m_nwindow_max) + ntable + " = " + string_utils::memsize(size) + "]";
}

std::string Buffer::capacity_string() const {
    return capacity_string(size());
}