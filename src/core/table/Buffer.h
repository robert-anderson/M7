//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include "src/core/util/defs.h"

class Buffer {
    template<typename table_t> friend
    class BufferedTable;

    template<typename table_t> friend
    class BufferedTableArray;

    std::vector <defs::data_t> m_data;

public:
    Buffer();

    Buffer(size_t dsize);

    Buffer(size_t ndword, size_t nrow);

    defs::data_t *ptr();

    const defs::data_t *ptr() const;

    size_t size() const;

    void resize(size_t dsize);

    void resize(size_t ndword, size_t nrow);
};



#endif //M7_BUFFER_H
