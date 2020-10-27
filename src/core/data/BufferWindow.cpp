//
// Created by RJA on 26/10/2020.
//

#include "BufferWindow.h"

BufferWindow::BufferWindow() : m_ptr(nullptr), m_dsize(0ul){}

BufferWindow::BufferWindow(Buffer &buffer) : m_ptr(buffer.ptr()), m_dsize(buffer.size()){}

BufferWindow::BufferWindow(Buffer &buffer, size_t doffset, size_t dsize) : m_ptr(buffer.ptr()+doffset), m_dsize(dsize){}

BufferWindow::operator bool() const {
    return m_ptr;
}
