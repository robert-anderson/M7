//
// Created by Robert J. Anderson on 26/10/2020.
//

#ifndef M7_BUFFER_H
#define M7_BUFFER_H

#include <list>

#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/parallel/SharedArray.h>

struct TableBase;

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

    /**
     * class representing the portion of a buffer allotted to a single table
     */
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

    private:
        /**
         * pointer to the beginning of the window
         */
        buf_t *m_begin_ptr = nullptr;
        /**
         * "high water mark" is a pointer to the highest-index in-use row
         */
        buf_t* m_hwm_ptr = nullptr;
        /**
         * while a window on a node-shared Buffer can be read (at the same pointer) across all MPI ranks on a node, it
         * must not be written to simultaneously by multiple MPI ranks on the same node because this would result in a
         * race condition. to prevent such situations, this "trash" buffer is provided for the non-root-on-node ranks to
         * write on, so that there is no need for the writing operations of Fields to be done within a conditional of
         * the form "if (mpi::on_node_i_am()) {...}"
         */
        v_t<buf_t> m_trash_buf;
    public:
        /**
         * number of bytes currently allotted to the window
         */
        uint_t m_size = 0ul;

        Window(uint_t row_size=1): m_row_size(row_size), m_trash_buf(m_row_size) {}

        Window(const Window& other): Window(other.m_row_size) {}

        Window& operator=(const Window& other);

        Window(Buffer *buffer, uint_t row_size=1);

        /**
         * @return
         *  pointer to the byte at the beginning of the buffer window with modifying intent
         */
        buf_t* begin() {return m_begin_ptr;}
        /**
         * @return
         *  pointer to the byte at the beginning of the buffer window with read-only intent
         */
        const buf_t* cbegin() const {return m_begin_ptr;}
        /**
         * @return
         *  pointer to the byte after the end of the buffer window with read-only intent
         */
        const buf_t* cend() const {return m_hwm_ptr;}

        /**
         * redefine the "high water mark"
         * @param irow
         *  index of the new "high water mark", must be no larger than the current number of allocated rows
         */
        void set_end(uint_t irow);

        /**
         * @return
         *  true if the buffer's underlying memory is shared over the MPI node
         */
        bool node_shared() const;

        /**
         * @return
         *  true if the memory is not node-shared, or if this is a node-root MPI rank
         */
        bool i_can_modify() const;

        /**
         * @return
         *  pointer to beginning of trash memory
         */
        buf_t* trash_dst() {return m_trash_buf.data();}

        bool operator==(const Window& other) const;

        /**
         * @return
         *  true if this window has an associated Buffer and it has an allocated data vector
         */
        bool allocated() const;

        /**
         * @return
         *  number of bytes in the usable range
         */
        uint_t size_in_use() const {
            return std::distance(cbegin(), cend());
        }

        bool empty() const {
            return cend() == cbegin();
        }

        void clear();

        /**
         * moves data if there's any in the window currently, and redefines the stored m_nbyte, m_nrow, and m_begin
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
        str_t name() const;

        void rename(str_t name) const {
            REQUIRE_TRUE(m_buffer, "only buffers can be named, this buffer window is not associated with one");
            m_buffer->m_name = std::move(name);
        }

        double get_expansion_factor() const;
    };


    double m_expansion_factor = 0.0;
    /**
     * if the buffer is given a name, its resizing events will be appear in the logs
     */
    str_t m_name = "";
    const uint_t m_nwindow_max;
    /**
     * determines whether the buffer is held in node-shared or rank-private memory
     */
    const bool m_node_shared;
private:
    /**
     * begin pointer of the allocated memory
     */
    buf_t* m_data = nullptr;
    /**
     * size of the memory pointed to by m_data
     */
    uint_t m_size = 0ul;
    /**
     * in the case that the buffer is initialized with m_shared=false
     */
    v_t<buf_t> m_data_priv;
    /**
     * in the case that the buffer is initialized with m_shared=true
     */
    SharedArrayBase m_data_shared;
    /**
     * a buffer can provide the underlying data requirement of multiple Tables, whose allocations are specified by
     * instances of the Window class
     */
    v_t<Window *> m_windows;

public:
    Buffer(str_t name, uint_t nwindow_max, bool node_shared=false);

    Buffer(uint_t nwindow_max, bool node_shared=false) : Buffer("", nwindow_max, node_shared){}

    uint_t size() const;

    uint_t window_size() const;

    void append_window(Window *window);

    /**
     * allocate a new data vector and call move method of all associated windows
     * @param new_size
     *  new new_size of the entire buffer in bytes
     * @param factor
     *  optional expansion factor
     */
    void resize(uint_t new_size, double factor=-1.0);
    /**
     * @param size
     *  size in bytes
     * @return
     *  space in memory taken up by the buffer in sensible units given the magnitude of the allocation
     */
    str_t capacity_string(uint_t size) const;

    str_t capacity_string() const;

};

struct Sizing {
    /**
     * estimate of the number of rows that will ultimately be required in a buffer
     */
    const uint_t m_nrec_est;
    /**
     * fractional over-allocation to make when buffer is resized
     */
    const double m_exp_fac;
};

#endif //M7_BUFFER_H
