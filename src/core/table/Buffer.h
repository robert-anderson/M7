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
        defs::data_t *m_dbegin = nullptr;
        defs::data_t *m_dend = nullptr;

        Window() {}

        explicit Window(Buffer *buffer);

        bool allocated() const;

        size_t dsize() const;

        void move(defs::data_t *dbegin, defs::data_t *dend);

        void resize(size_t dsize);

        void make_room(size_t dsize);

        void expand(size_t delta_dsize, double expansion_factor);

        void expand(size_t delta_dsize);

        double expansion_factor() const;

        std::string name() const {
            return m_buffer->m_name;
        }

    };

    double m_expansion_factor = 0.5;
    const std::string m_name;
private:
    const size_t m_nwindow_max;
    std::vector<defs::data_t> m_data;
    std::vector<Window *> m_windows;


public:
    Buffer(std::string name, size_t nwindow_max);

    size_t dsize() const;

    size_t window_dsize() const;

    void append_window(Window *window);

    void resize(size_t dsize);

    // resize if smaller
    void make_room(size_t dsize);

    /**
     *
     * @param delta_nrow
     * number of rows to add to the buffer
     * @param expansion_factor
     * must be non-negative. If zero, the buffer is expanded
     * by delta_nrow rows exactly. Any other value indicates the proportion of
     * additional rows to add relative to the new total number of rows.
     */
    void expand(size_t delta_dsize, double expansion_factor);

    void expand(size_t delta_dsize);

    std::string capacity_string(size_t dsize) const;

    std::string capacity_string() const;

};

#endif //M7_BUFFER_H
