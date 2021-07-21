//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include <list>
#include "src/defs.h"
#include "src/core/parallel/MPIWrapper.h"


class Buffer {

public:
    class Window {
        friend class Buffer;

        Buffer *m_buffer = nullptr;

    public:
        defs::buf_t *m_begin = nullptr;
        defs::buf_t *m_end = nullptr;

        Window() {}

        explicit Window(Buffer *buffer);

        bool allocated() const;

        size_t size() const;

        void move(defs::buf_t *begin, defs::buf_t *end);

        void resize(size_t size, double factor=-1.0);

        std::string name() const;

        double get_expansion_factor() const;

    };

    double m_expansion_factor = 0.0;
    const std::string m_name;
private:
    const size_t m_nwindow_max;
    std::vector<defs::buf_t> m_data;
    std::vector<Window *> m_windows;

public:
    Buffer(std::string name, size_t nwindow_max);

    size_t size() const;

    size_t window_size() const;

    void append_window(Window *window);

    /**
     * allocate a
     * @param size
     * @param factor
     */
    void resize(size_t size, double factor=-1.0);

    std::string capacity_string(size_t dsize) const;

    std::string capacity_string() const;

};

#endif //M7_BUFFER_H
