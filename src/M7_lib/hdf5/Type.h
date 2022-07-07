//
// Created by Robert J. Anderson on 13/12/2020.
//

#ifndef M7_HDF5_TYPE_H
#define M7_HDF5_TYPE_H

#include <string>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/nd/NdFormat.h>

namespace hdf5 {
#ifdef H5_HAVE_PARALLEL
    constexpr bool c_have_parallel = true;
#else
    constexpr bool c_have_parallel = false;
#endif

    static_assert(c_have_parallel, "HDF5 must be compiled with parallel functionality");



    static const std::array<hid_t, 14> c_types = {
        0, H5T_NATIVE_CHAR, H5T_NATIVE_SCHAR, H5T_NATIVE_UCHAR, H5T_NATIVE_SHORT, H5T_NATIVE_USHORT,
        H5T_NATIVE_INT, H5T_NATIVE_UINT, H5T_NATIVE_LONG, H5T_NATIVE_ULONG, H5T_NATIVE_LLONG, H5T_NATIVE_ULLONG,
        H5T_NATIVE_FLOAT, H5T_NATIVE_DOUBLE
    };

    str_t type_name(hid_t type);

    template<typename T=void>
    static constexpr uint_t type_ind() { return 0; }

    template<> constexpr uint_t type_ind<char>() { return 1; }
    template<> constexpr uint_t type_ind<signed char>() { return 2; }
    template<> constexpr uint_t type_ind<unsigned char>() { return 3; }
    template<> constexpr uint_t type_ind<short>() { return 4; }
    template<> constexpr uint_t type_ind<unsigned short>() { return 5; }
    template<> constexpr uint_t type_ind<int>() { return 6; }
    template<> constexpr uint_t type_ind<unsigned int>() { return 7; }
    template<> constexpr uint_t type_ind<long>() { return 8; }
    template<> constexpr uint_t type_ind<unsigned long>() { return 9; }
    template<> constexpr uint_t type_ind<long long>() { return 10; }
    template<> constexpr uint_t type_ind<unsigned long long>() { return 11; }
    template<> constexpr uint_t type_ind<float>() { return 12; }
    template<> constexpr uint_t type_ind<double>() { return 13; }


    /**
     * @param h5type
     *  type index
     * @return
     *  size in bytes of the type identified by the given HDF5 type index
     */
    hsize_t type_size(hid_t h5type);



    struct Type {
        const hid_t m_handle;
        const hsize_t m_size;

    private:
        const bool m_immutable;
        static hsize_t size_max(const strv_t* vec);
        /*
         * use dummy arg so as not to have same prototype as public ctor in case hsize_t coincides with hid_t
         */
        Type(hsize_t size, char /*dummy*/);

    public:

        template<typename T>
        static const hid_t &type() {
            typedef arith::comp_t<T> comp_t;
            static_assert(type_ind<comp_t>(), "type has no HDF5 native equivalent");
            return c_types[type_ind<comp_t>()];
        }

        Type(): m_handle(0), m_size(0ul), m_immutable(true){}
        explicit Type(hid_t handle): m_handle(handle), m_size(H5Tget_size(m_handle)), m_immutable(true){}

        template<typename T>
        Type(const T*): Type(type<T>()){}


        Type(const str_t* str);

        Type(const strv_t* str_vec);

        ~Type();

        static bool types_equal(hid_t t1, hid_t t2) {
            return H5Tequal(t1, t2);
        }

        bool operator==(hid_t other) const {
            return types_equal(m_handle, other);
        }

        operator hid_t () const;

    };
}


#endif //M7_HDF5_TYPE_H
