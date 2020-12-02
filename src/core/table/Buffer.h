//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include "src/defs.h"

class Buffer {
    template<typename table_t> friend
    class BufferedTable;

    template<typename table_t> friend
    class BufferedTableArray;

    std::string m_name;
    size_t m_row_dsize, m_nrow;
    std::vector<defs::data_t> m_data;

public:

    Buffer(std::string name, size_t row_dsize, size_t nrow);

    ~Buffer();

    Buffer(const Buffer &other);

    Buffer(Buffer &&other);

    Buffer &operator=(const Buffer &other);

    Buffer &operator=(Buffer &&other);

    defs::data_t *ptr();

    const defs::data_t *ptr() const;

    size_t dsize() const;

    std::string freeing_string() const;

    std::string capacity_string() const;

    void resize(size_t row_dsize, size_t nrow);
};


#endif //M7_BUFFER_H
