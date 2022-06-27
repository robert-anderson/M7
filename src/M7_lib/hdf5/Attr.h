//
// Created by anderson on 27/06/2022.
//

#ifndef M7_HDF5_ATTR_H
#define M7_HDF5_ATTR_H

#include "Dataspace.h"

namespace hdf5 {
    class AttrReader {
        const hid_t m_handle;
        const DataSpace m_space;
    public:
        const Type m_type;
        const hsize_t m_nelement;

        AttrReader(hid_t parent_handle, const std::string &name);

        ~AttrReader();

        void read_bytes(char *dst) const;

        template<typename T>
        void read(T *dst) const {
            REQUIRE_EQ(Type(dst), m_type, "element type is at odds with the stored type");
            read_bytes(reinterpret_cast<char *>(dst));
        }

        void read(std::string *dst, size_t n) const;
    };

    class AttrWriter {
        const DataSpace m_space;
        const hid_t m_h5type;
        const hid_t m_handle;

    public:
        AttrWriter(hid_t parent_handle, const std::string &name, const std::vector <hsize_t> &shape, hid_t h5type);

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
