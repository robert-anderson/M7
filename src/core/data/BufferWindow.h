//
// Created by rja on 02/10/2020.
//

#ifndef M7_BUFFERWINDOW_H
#define M7_BUFFERWINDOW_H


#include "Buffer.h"

struct BufferWindow {
    // pointer to the initial data_t element
    defs::data_t* m_ptr;
    // number of data_t elements allotted to this window
    size_t m_dsize;

    BufferWindow();

    BufferWindow(Buffer& buffer);

    BufferWindow(Buffer &buffer, size_t doffset, size_t dsize);

    operator bool() const;
};



#endif //M7_BUFFERWINDOW_H
