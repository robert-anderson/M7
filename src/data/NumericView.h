//
// Created by Robert John Anderson on 2020-02-12.
//

#ifndef M7_NUMERICVIEW_H
#define M7_NUMERICVIEW_H

#include <src/utils.h>

template<typename T>
class NumericView {
    T *const m_data;
    const size_t m_size;
public:
    NumericView(const char *data, const size_t size = 1) :
            m_data((T *) data), m_size(size) {}

    void operator=(const T &value) { m_data[0] = value; }

    T &operator[](size_t i) const {
        assert(i < m_size);
        return m_data[i];
    }

    T &operator*() const { return *m_data; }

    const size_t size() const { return m_size; }

    std::string to_string(size_t padding = 0) const {
        std::string out{};
        for (auto i{0ul}; i < size(); ++i) {
            out += utils::num_to_string(m_data[i], padding);
        }
        return out;
    }
};


#endif //M7_NUMERICVIEW_H
