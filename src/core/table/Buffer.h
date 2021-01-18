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
        defs::data_t *m_dbegin = nullptr;
        defs::data_t *m_dend = nullptr;

    public:
        Window() {}

        explicit Window(Buffer *buffer);

        size_t dsize() const;

        size_t size() const;

        void move(defs::data_t *dbegin, defs::data_t *dend);

        defs::data_t *dbegin();

        const defs::data_t *dbegin() const;

        void resize(size_t dsize);

        void expand();

        void expand(size_t delta_dsize);

        double expansion_factor() const;

    };

    double m_expansion_factor = 0.5;
    const std::string m_name;
private:
    const size_t m_nwindow_max;
    std::vector<defs::data_t> m_data;
    std::vector<Window *> m_windows;

public:
    Buffer(std::string name, size_t nwindow_max);

    Buffer(std::string name, size_t nwindow_max, size_t dsize);

    size_t dsize() const;

    size_t window_dsize() const;

    defs::data_t *dbegin();

    const defs::data_t *dbegin() const;

    defs::data_t *dbegin(const size_t &iwindow);

    const defs::data_t *dbegin(const size_t &iwindow) const;

    void append_window(Window *window);

    void resize(size_t dsize);

    void expand();

    void expand(size_t delta_dsize);;

    std::string capacity_string() const;

};

#endif //M7_BUFFER_H
