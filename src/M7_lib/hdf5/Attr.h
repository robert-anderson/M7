//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_ATTR_H
#define M7_HDF5_ATTR_H

#include <utility>

#include "Dataspace.h"
#include "IoManager.h"

namespace hdf5 {

    struct Attr {
        const v_t<buf_t> m_buf;
        const dataset::Format m_format;
        Attr(v_t<buf_t> buf, dataset::Format format): m_buf(std::move(buf)), m_format(std::move(format)){
            REQUIRE_EQ(m_buf.size(), m_format.m_size, "buffer size inconsistent with format");
        }

        template<typename T>
        Attr(const T* v, uint_t size):
            Attr({reinterpret_cast<const buf_t*>(v), reinterpret_cast<const buf_t*>(v + size)},
                 dataset::Format(Type::make<T>(), {size}, {}, dtype::is_complex<T>())){}

        template<typename T>
        Attr(const T& v): Attr(&v, 1ul){}

        template<typename T>
        Attr(const v_t<T>& v): Attr(v.data(), v.size()){}

    private:

        template<typename T>
        bool parsable_as() {
            return Type::make<T>() == m_format.m_h5_type;
        }

        template<typename T>
        bool parse(T* v) {
            if (!parsable_as<T>()) return false;
            auto dst = reinterpret_cast<buf_t*>(v);
            std::memcpy(dst, m_buf.data(), m_buf.size());
            return true;
        }

    public:
        template<typename T>
        bool parse(v_t<T>& v) {
            if (!parsable_as<T>()) return false;
            DEBUG_ASSERT_FALSE(m_buf.size() % sizeof(T), "buffer size is not integer multiple of parsing type size");
            v.resize(m_buf.size() / sizeof(T));
            return parse(v.data());
        }

        template<typename T>
        bool parse(T& v) {
            REQUIRE_EQ(sizeof(T), m_buf.size(), "buffer size is equal to that of a single word of the parsing type");
            return parse(&v);
        }
    };

    class AttrReader {
        const hid_t m_handle;
        const DataSpace m_space;
    public:
        const Type m_type;
        const hsize_t m_nelement;

        AttrReader(hid_t parent_handle, const str_t &name);

        ~AttrReader();

        void read_bytes(char *dst) const;

        template<typename T>
        void read(T *dst) const {
            REQUIRE_EQ(Type(dst), m_type, "element type is at odds with the stored type");
            read_bytes(reinterpret_cast<char *>(dst));
        }

        void read(str_t *dst, size_t n) const;
    };

    class AttrWriter {
        const DataSpace m_space;
        const hid_t m_h5type;
        const hid_t m_handle;

    public:
        AttrWriter(hid_t parent_handle, const str_t &name, const v_t <hsize_t> &shape, hid_t h5type);

        ~AttrWriter();

        void write_bytes(const char *src) const;

        template<typename T>
        void write(const T *src, hsize_t n) const {
            REQUIRE_EQ(Type(src), m_h5type, "element type is at odds with the stored type");
            REQUIRE_EQ(n, m_space.m_nelement, "number of elements written must match that of the dataspace");
            write_bytes(reinterpret_cast<const char *>(src));
        }
    };
}

#endif //M7_HDF5_ATTR_H
