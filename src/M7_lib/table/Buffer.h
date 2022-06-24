//
// Created by Robert J. Anderson on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include <list>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/parallel/MPIWrapper.h>


class Buffer {
public:
    /**
     * number of bytes in the system word
     */
    static constexpr uint_t c_nbyte_word = sizeof(uint_t);
    /**
     * number of bits in the system word
     */
    static constexpr uint_t c_nbit_word = CHAR_BIT * c_nbyte_word;
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
        const uint_t m_row_size;
        /**
         * Current number of whole rows that can be stored in the window
         */
        uint_t m_nrow = 0ul;
        buf_t *m_begin = nullptr;
        uint_t m_size = 0ul;

        Window(uint_t row_size=1): m_row_size(row_size) {}

        Window(const Window& other): Window(other.m_row_size) {}

        Window& operator=(const Window& other);

        explicit Window(Buffer *buffer, uint_t row_size=1);
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
        void move(buf_t *begin, uint_t new_size);
        /**
         * delegate to the buffer's resize method so that all associated windows, not just this one, are moved
         * @param size
         *  minimum required new size in bytes
         * @param factor
         *  optional expansion factor
         */
        void resize(uint_t size, double factor=-1.0);
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
    const uint_t m_nwindow_max;
    std::vector<buf_t> m_data;
    std::vector<Window *> m_windows;

public:
    Buffer(std::string name, uint_t nwindow_max);

    uint_t size() const;

    uint_t window_size() const;

    void append_window(Window *window);

    /**
     * allocate a new data vector and call move method of all associated windows
     * @param size
     *  new size of the entire buffer in bytes
     * @param factor
     *  optional expansion factor
     */
    void resize(uint_t size, double factor=-1.0);
    /**
     * @param size
     *  size in bytes
     * @return
     *  space in memory taken up by the buffer in sensible units given the magnitude of the allocation
     */
    std::string capacity_string(uint_t size) const;

    std::string capacity_string() const;

};

#endif //M7_BUFFER_H
