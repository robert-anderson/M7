//
// Created by Robert John Anderson on 2020-02-12.
//

#ifndef M7_NUMERICVIEW_H
#define M7_NUMERICVIEW_H

#include <src/utils.h>
#include <cstring>

template<typename T>
class NumericView {
    T *const m_data;
    const size_t m_size;
public:
    NumericView(const char *data, const size_t size = 1) :
            m_data((T *) data), m_size(size) {}

    void operator=(const T &value) { m_data[0] = value; }
    void operator=(const NumericView<T> &value) { m_data[0] = value.m_data[0]; }

    T &operator[](size_t i) const {
        assert(i < m_size);
        return m_data[i];
    }

    T &operator*() const { return *m_data; }

    void operator*=(const T& value){m_data[0]*=value;}
    void operator+=(const T& value){m_data[0]+=value;}
    void operator+=(const NumericView<T>& value){m_data[0]+=value.m_data[0];}


    // allow the view to be implicitly converted to the data type it holds
    operator T() const {
        assert(m_data!= nullptr);
        return *m_data;
    }

    const size_t size() const { return m_size; }

    std::string to_string(size_t padding = 0) const {
        std::string out{};
        if (!m_data) return out;
        for (size_t i=0ul; i < size(); ++i) {
            out += utils::num_to_string(m_data[i], padding);
        }
        return out;
    }

    bool is_zero() {
        return this == 0;
    }

    void zero() {
        memset((void *) m_data, 0, m_size * sizeof(T));
    }
};


#endif //M7_NUMERICVIEW_H
