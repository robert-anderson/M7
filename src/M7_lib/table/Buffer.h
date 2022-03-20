//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include <list>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/parallel/MPIWrapper.h>


class Buffer {

public:
    class Window {
        friend class Buffer;
        /**
         * pointer to the Buffer with which this window is associated
         */
        Buffer *m_buffer = nullptr;

    public:
        /**
         * Size of the row in bytes
         */
        const size_t m_row_size;
        /**
         * Current number of whole rows that can be stored in the window
         */
        size_t m_nrow = 0ul;
        defs::buf_t *m_begin = nullptr;
        size_t m_size = 0ul;

        Window(size_t row_size=1): m_row_size(row_size) {}

        Window(const Window& other): Window(other.m_row_size) {}

        Window& operator=(const Window& other);

        explicit Window(Buffer *buffer, size_t row_size=1);
        /**
         * @return
         *  true if this window has an associated Buffer and it has an allocated data vector
         */
        bool allocated() const;
        /**
         * moves data if there's any in the window currently, and redefines the stored m_size, m_nrow, and m_begin
         * @param begin
         *  new begin buffer pointer
         * @param new_size
         *  new size in bytes for the redefined window
         */
        void move(defs::buf_t *begin, size_t new_size);
        /**
         * delegate to the buffer's resize method so that all associated windows, not just this one, are moved
         * @param size
         *  minimum required new size in bytes
         * @param factor
         *  optional expansion factor
         */
        void resize(size_t size, double factor=-1.0);
        /**
         * @return
         *  name of buffer
         */
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
     * allocate a new data vector and call move method of all associated windows
     * @param size
     *  new size of the entire buffer in bytes
     * @param factor
     *  optional expansion factor
     */
    void resize(size_t size, double factor=-1.0);
    /**
     * @param size
     *  size in bytes
     * @return
     *  space in memory taken up by the buffer in sensible units given the magnitude of the allocation
     */
    std::string capacity_string(size_t size) const;

    std::string capacity_string() const;

};

#endif //M7_BUFFER_H
