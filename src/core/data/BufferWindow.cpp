//
// Created by rja on 02/10/2020.
//

#include "BufferWindow.h"

BufferWindow::BufferWindow() : m_ptr(nullptr), m_dsize(0ul){}

BufferWindow::BufferWindow(Buffer &buffer) : m_ptr(buffer.ptr()), m_dsize(buffer.size()){}

BufferWindow::operator bool() const {
    return m_ptr;
}
