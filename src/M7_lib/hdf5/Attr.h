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
        const dataset::ItemFormat m_format;
        const str_t m_name;
    private:
        /*
         * char arg is to resolve ambiguity between this ctor definition and the next in the case that T=buf_t
         */
        Attr(v_t<buf_t> buf, dataset::ItemFormat format, str_t name, char /*dummy*/);
    public:

        template<typename T>
        Attr(const T* buf, dataset::ItemFormat format, str_t name):
            Attr({reinterpret_cast<const buf_t*>(buf), reinterpret_cast<const buf_t*>(buf)+format.m_size}, format, name, 0){}

        template<typename T>
        Attr(const v_t<T>& buf, const uintv_t& shape, str_t name):
                Attr(buf.data(), {Type::make<T>(), shape, {}, dtype::is_complex<T>()}, name){
            REQUIRE_EQ(buf.size(), nd::nelement(shape), "format is not compatible with given vector");
        }

        Attr(hid_t parent_handle, str_t name);

        template<typename T>
        Attr(const T* v, uint_t size, str_t name): Attr(v, {Type(v), {size}, {}, dtype::is_complex<T>()}, name){}

        /*
         * do not include the null terminator in the length of the stored string
         */
        Attr(const str_t& v, str_t name): m_buf(v.data(), v.data()+v.size()),
            m_format(Type(&v), {1}, {}, false), m_name(std::move(name)){}

        template<typename T>
        Attr(const T& v, str_t name): Attr(&v, 1ul, std::move(name)){}

        template<typename T>
        Attr(const v_t<T>& v, str_t name): Attr(v.data(), v.size(), std::move(name)){}

        template<typename T, uint_t nind>
        Attr(const std::array<T, nind>& a, str_t name): Attr(convert::to_vector(a.data(), nind), std::move(name)){}

        bool operator==(const Attr& other) const;

        bool operator!=(const Attr& other) const;

    private:

        template<typename T>
        bool parsable_as() {
            if (dtype::is_complex<T>() && (m_format.m_h5_shape.back()!=2ul))
                // minor dimension must be 2, for the real/imag components
                return false;
            return Type::make<T>() == m_format.m_type;
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

        template<typename T>
        void parse(T& v, const T& default_) {
            if (!parsable_as<T>()) v = default_;
            else parse(v);
        }

        void save(hid_t parent_handle) const;
        static Attr load(hid_t parent_handle, const str_t& name);
    };

}

#endif //M7_HDF5_ATTR_H
